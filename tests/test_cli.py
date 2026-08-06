from gpcr_selectivitymap.cli import parser


def test_cli_commands_exist() -> None:
    p = parser()
    assert p.prog == "gpcr-selectivitymap"
    choices = p._subparsers._group_actions[0].choices
    for command in (
        "profile",
        "batch",
        "benchmark",
        "cluster-benchmark",
        "phylogeny-baseline",
        "interface-audit",
        "provenance-audit",
        "aggregate",
        "complementarity-audit",
        "receptor-compatibility",
    ):
        assert command in choices
