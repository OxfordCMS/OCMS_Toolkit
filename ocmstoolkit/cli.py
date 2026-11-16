"""Command-line interface for OCMS_Toolkit.

Each Snakemake workflow shipped in ``ocmstoolkit/workflows/pipeline_<name>.smk``
is exposed as a sub-command ``<name>``, with a matching configuration template
in ``ocmstoolkit/config/pipeline_<name>.yaml``.

Typical use, from the directory in which the pipeline should run::

    ocms_toolkit fastqc config            # write ./pipeline_fastqc.yaml
    # ... edit pipeline_fastqc.yaml ...
    ocms_toolkit fastqc full --dry        # show the jobs that would run
    ocms_toolkit fastqc full --jobs 20    # submit the jobs to Slurm

Options after a bare ``--`` are passed unchanged to Snakemake, for example::

    ocms_toolkit fastqc full -- --unlock

The scripts in ``ocmstoolkit/scripts`` are run with their own options::

    ocms_toolkit rename_and_link --help
"""

from __future__ import annotations

import argparse
import getpass
import importlib
import importlib.metadata
import importlib.util
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

PROG = "ocms_toolkit"

PACKAGE_DIR = Path(__file__).resolve().parent
WORKFLOW_DIR = PACKAGE_DIR / "workflows"
CONFIG_DIR = PACKAGE_DIR / "config"
SCRIPT_DIR = PACKAGE_DIR / "scripts"

PIPELINE_PREFIX = "pipeline_"

#: Value used in the configuration templates until a real account is set.
ACCOUNT_PLACEHOLDER = "your_slurm_account"

#: One-line descriptions of the scripts, shown in ``ocms_toolkit --help``.
SCRIPT_SUMMARIES = {
    "fastq2filteredfastq": "filter a FASTQ file by the reads mapped (or "
                           "unmapped) in a BAM file",
    "new_project": "set up the directories for a new project on BMRC",
    "rename_and_link": "symlink files under new names from a mapping file",
}


class CLIError(Exception):
    """An error reported to the user without a traceback."""


# ---------------------------------------------------------------------------
# Pipeline and script discovery
# ---------------------------------------------------------------------------


def discover_pipelines(workflow_dir: Path = WORKFLOW_DIR) -> dict[str, Path]:
    """Return ``{pipeline name: Snakefile}`` for every ``pipeline_*.smk``."""
    return {
        smk.stem[len(PIPELINE_PREFIX):]: smk
        for smk in sorted(workflow_dir.glob(f"{PIPELINE_PREFIX}*.smk"))
    }


def discover_scripts(script_dir: Path = SCRIPT_DIR) -> list[str]:
    """Names of the scripts in ``ocmstoolkit/scripts``."""
    return sorted(
        path.stem for path in script_dir.glob("*.py")
        if not path.stem.startswith("_")
    )


def config_filename(name: str) -> str:
    """File name of the configuration for pipeline ``name``."""
    return f"{PIPELINE_PREFIX}{name}.yaml"


def config_template(name: str, config_dir: Path = CONFIG_DIR) -> Path:
    """Packaged configuration template for pipeline ``name``."""
    return config_dir / config_filename(name)


def pipeline_summary(snakefile: Path) -> str:
    """First comment line of a Snakefile, used as the sub-command help."""
    try:
        with snakefile.open() as handle:
            for line in handle:
                line = line.strip()
                if line.startswith("#"):
                    return line.lstrip("#").strip()
                if line:
                    break
    except OSError:
        pass
    return f"Run the {snakefile.stem[len(PIPELINE_PREFIX):]} pipeline"


# ---------------------------------------------------------------------------
# Slurm account handling
# ---------------------------------------------------------------------------


def _sacctmgr(*args: str) -> list[str] | None:
    """Run ``sacctmgr`` and return its output lines.

    Returns ``None`` when ``sacctmgr`` is unavailable or fails, so that callers
    can distinguish "cannot check" from "no results".
    """
    exe = shutil.which("sacctmgr")
    if exe is None:
        return None
    try:
        result = subprocess.run(
            [exe, "--noheader", "--parsable2", *args],
            capture_output=True,
            text=True,
            check=True,
            timeout=60,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def default_slurm_account(user: str) -> str | None:
    """The user's default Slurm account, or ``None`` if it cannot be found."""
    rows = _sacctmgr("show", "user", user, "format=DefaultAccount")
    return rows[0] if rows else None


def slurm_accounts(user: str) -> set[str] | None:
    """All Slurm accounts the user may submit to, or ``None`` if unknown."""
    rows = _sacctmgr("show", "associations", f"user={user}", "format=Account")
    return None if rows is None else set(rows)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


def write_config(
    name: str,
    destination: Path,
    force: bool = False,
    account: str | None = None,
    config_dir: Path = CONFIG_DIR,
) -> bool:
    """Copy the configuration template for ``name`` to ``destination``.

    If ``account`` is given it replaces the Slurm account placeholder. The
    template is copied as text so that its comments are preserved.

    Returns ``True`` if the account placeholder was filled in.
    """
    template = config_template(name, config_dir)
    if not template.is_file():
        raise CLIError(
            f"No configuration template for pipeline '{name}' "
            f"(expected {template})."
        )
    if destination.exists() and not force:
        raise CLIError(
            f"{destination} already exists. Edit it, or use --force to "
            "replace it with a fresh copy of the template."
        )

    text = template.read_text()
    filled = False
    if account:
        text, count = re.subn(
            rf"^(\s*slurm_account:\s*)[\"']?{re.escape(ACCOUNT_PLACEHOLDER)}[\"']?",
            lambda m: f'{m.group(1)}"{account}"',
            text,
            flags=re.MULTILINE,
        )
        filled = count > 0
    destination.write_text(text)
    return filled


def load_config(path: Path) -> dict:
    """Read a pipeline configuration file."""
    try:
        with path.open() as handle:
            config = yaml.safe_load(handle)
    except yaml.YAMLError as err:
        raise CLIError(f"{path} is not valid YAML:\n{err}") from err
    if not isinstance(config, dict):
        raise CLIError(f"{path} does not contain a YAML mapping.")
    return config


def get_resources(config: dict) -> dict:
    """The ``resources`` section of a configuration (may be empty)."""
    resources = config.get("resources") or {}
    if not isinstance(resources, dict):
        raise CLIError("'resources' in the configuration must be a mapping.")
    return resources


def validate_slurm(config: dict, user: str) -> list[str]:
    """Check the Slurm settings needed to submit jobs.

    Raises :class:`CLIError` for problems that would make submission fail, and
    returns a list of warnings for checks that could not be carried out.
    """
    resources = get_resources(config)
    account = resources.get("slurm_account")
    if not account or str(account) == ACCOUNT_PLACEHOLDER:
        raise CLIError(
            "No Slurm account is set. Set resources.slurm_account in the "
            "configuration file (list your accounts with: "
            f"sacctmgr show associations user={user} format=Account)."
        )

    warnings = []
    accounts = slurm_accounts(user)
    if accounts is None:
        warnings.append(
            "Could not query sacctmgr; the Slurm account "
            f"'{account}' has not been verified."
        )
    elif str(account) not in accounts:
        available = ", ".join(sorted(accounts)) or "none found"
        raise CLIError(
            f"'{account}' is not one of your Slurm accounts "
            f"(available: {available})."
        )
    return warnings


def format_default_resources(resources: dict) -> list[str]:
    """Render ``resources`` as ``--default-resources`` arguments.

    Snakemake evaluates each value as a Python expression, so text values are
    quoted and numbers are passed unchanged.
    """
    rendered = []
    for key, value in resources.items():
        if value is None or value == "":
            continue
        if isinstance(value, bool) or not isinstance(value, (int, float, str)):
            raise CLIError(
                f"resources.{key} must be a number or text, not {value!r}."
            )
        if isinstance(value, str):
            rendered.append(f"{key}={value!r}")
        else:
            rendered.append(f"{key}={value}")
    return rendered


# ---------------------------------------------------------------------------
# Snakemake invocation
# ---------------------------------------------------------------------------


def build_snakemake_command(
    snakefile: Path,
    configfile: Path,
    config: dict,
    jobs: int,
    cores: int | None = None,
    dry_run: bool = False,
    local: bool = False,
    extra: list[str] | None = None,
) -> list[str]:
    """Assemble the Snakemake command line for a pipeline run."""
    cmd = [
        sys.executable, "-m", "snakemake",
        "--snakefile", str(snakefile),
        "--configfile", str(configfile),
        "--directory", str(configfile.parent),
    ]
    if local:
        cmd += ["--cores", str(cores) if cores else "all"]
    else:
        # --cores is deliberately not passed to Slurm runs: Snakemake caps
        # the threads of every job at --cores, so "--cores 1" would submit
        # all jobs with --cpus-per-task=1. --local-cores only limits the
        # rules that run on this node.
        cmd += ["--executor", "slurm", "--jobs", str(jobs)]
        if cores:
            cmd += ["--local-cores", str(cores)]

    # Applied to every job, including rules that define their own
    # resources, so that the Slurm account always reaches sbatch.
    default_resources = format_default_resources(get_resources(config))
    if default_resources:
        cmd += ["--default-resources", *default_resources]

    cmd += ["--printshellcmds", "--rerun-incomplete"]
    if dry_run:
        cmd.append("--dry-run")
    cmd += list(extra or [])
    return cmd


def run_pipeline(name: str, snakefile: Path, args: argparse.Namespace,
                 extra: list[str]) -> int:
    """Validate the configuration and run Snakemake for ``name``."""
    if importlib.util.find_spec("snakemake") is None:
        raise CLIError(
            "Snakemake is not installed in this Python environment "
            f"({sys.executable})."
        )

    configfile = (Path.cwd() / config_filename(name)).resolve()
    if not configfile.is_file():
        raise CLIError(
            f"{configfile.name} not found in {configfile.parent}.\n"
            f"Create it with: {PROG} {name} config"
        )
    config = load_config(configfile)

    if not args.local:
        for warning in validate_slurm(config, getpass.getuser()):
            print(f"Warning: {warning}", file=sys.stderr)

    cmd = build_snakemake_command(
        snakefile=snakefile,
        configfile=configfile,
        config=config,
        jobs=args.jobs,
        cores=args.cores,
        dry_run=args.dry,
        local=args.local,
        extra=extra,
    )
    print("Running:", shlex.join(cmd), file=sys.stderr, flush=True)
    try:
        return subprocess.run(cmd).returncode
    except KeyboardInterrupt:
        return 130


def make_config(name: str, args: argparse.Namespace) -> int:
    """Write the configuration file for ``name`` to the current directory."""
    destination = Path.cwd() / config_filename(name)
    account = default_slurm_account(getpass.getuser())
    filled = write_config(name, destination, force=args.force, account=account)

    print(f"Wrote {destination}")
    if filled:
        print(f"Slurm account set to your default account '{account}'.")
    else:
        print("Set resources.slurm_account before running the pipeline.")
    print(
        f"Edit the file, then check the jobs with: {PROG} {name} full --dry"
    )
    return 0


def run_script(name: str, argv: list[str]) -> int:
    """Run ``ocmstoolkit/scripts/<name>.py`` with the arguments ``argv``."""
    module = importlib.import_module(f"ocmstoolkit.scripts.{name}")
    script_argv = [f"{PROG} {name}", *argv]
    # Scripts that parse sys.argv themselves see the same arguments.
    saved_argv = sys.argv
    sys.argv = script_argv
    try:
        result = module.main(script_argv)
    finally:
        sys.argv = saved_argv
    return result if isinstance(result, int) else 0


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


def _version() -> str:
    try:
        return importlib.metadata.version("ocms_toolkit")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _scripts_epilog(scripts: list[str]) -> str:
    width = max((len(name) for name in scripts), default=0)
    lines = ["scripts (run '%(prog)s <script> --help' for their options):"]
    for name in scripts:
        summary = SCRIPT_SUMMARIES.get(name, "")
        lines.append(f"  {name:<{width}}  {summary}".rstrip())
    return "\n".join(lines)


def build_parser(pipelines: dict[str, Path],
                 scripts: list[str] | None = None) -> argparse.ArgumentParser:
    """Build the argument parser with one sub-command per pipeline."""
    epilog = (
        "Run '%(prog)s <pipeline> <task> --help' for the options of a task.\n"
        "Arguments after a bare '--' are passed to Snakemake."
    )
    if scripts:
        epilog += "\n\n" + _scripts_epilog(scripts)
    parser = argparse.ArgumentParser(
        prog=PROG,
        description="Run the OCMS_Toolkit Snakemake pipelines and scripts.",
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {_version()}")
    subparsers = parser.add_subparsers(dest="pipeline", metavar="<pipeline>",
                                       title="pipelines")

    for name, snakefile in pipelines.items():
        summary = pipeline_summary(snakefile)
        pipeline = subparsers.add_parser(name, help=summary, description=summary)
        tasks = pipeline.add_subparsers(dest="task", metavar="<task>",
                                        title="tasks", required=True)

        config = tasks.add_parser(
            "config",
            help=f"write {config_filename(name)} to the current directory",
            description=(
                f"Write {config_filename(name)} to the current directory, with "
                "your default Slurm account filled in where it can be found."
            ),
        )
        config.add_argument("--force", action="store_true",
                            help="replace an existing configuration file")

        full = tasks.add_parser(
            "full",
            help="run the complete pipeline",
            description=(
                f"Run the complete {name} pipeline using "
                f"{config_filename(name)} in the current directory."
            ),
        )
        full.add_argument("-j", "--jobs", type=int, default=10,
                          help="maximum number of concurrent Slurm jobs "
                               "(default: %(default)s)")
        full.add_argument("-c", "--cores", type=int, default=None,
                          help="with --local, the number of cores to use "
                               "(default: all); otherwise the cores for the "
                               "steps that run on this node (default: all)")
        full.add_argument("-n", "--dry", action="store_true",
                          help="show the jobs that would run, without running")
        full.add_argument("--local", action="store_true",
                          help="run on the current node instead of "
                               "submitting to Slurm")
    return parser


def _cgat_syntax_hint(argv: list[str], pipelines: dict[str, Path]) -> str | None:
    """Explain the new syntax if the old cgat ``make`` syntax was used."""
    if len(argv) >= 2 and argv[0] in pipelines and argv[1] in ("make", "show"):
        return (
            f"'{PROG} {argv[0]} {argv[1]} ...' is the syntax of the earlier "
            "cgatcore pipelines. Use:\n"
            f"  {PROG} {argv[0]} full --dry       (show the jobs)\n"
            f"  {PROG} {argv[0]} full --jobs N    (run the pipeline)"
        )
    return None


def main(argv: list[str] | None = None) -> int:
    """Entry point for the ``ocms_toolkit`` command."""
    argv = list(sys.argv[1:] if argv is None else argv)

    scripts = discover_scripts()
    if argv and argv[0] in scripts:
        return run_script(argv[0], argv[1:])

    extra: list[str] = []
    if "--" in argv:
        split = argv.index("--")
        argv, extra = argv[:split], argv[split + 1:]

    pipelines = discover_pipelines()
    hint = _cgat_syntax_hint(argv, pipelines)
    if hint:
        print(f"Error: {hint}", file=sys.stderr)
        return 2

    parser = build_parser(pipelines, scripts)
    args = parser.parse_args(argv)

    if args.pipeline is None:
        parser.print_help()
        return 1
    if extra and args.task != "full":
        parser.error("arguments after '--' are only accepted by 'full'")

    try:
        if args.task == "config":
            return make_config(args.pipeline, args)
        return run_pipeline(args.pipeline, pipelines[args.pipeline], args, extra)
    except CLIError as err:
        print(f"Error: {err}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
