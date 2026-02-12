"""CLI entry point para ejecutar Chemuson como módulo."""


def main() -> None:
    from chemuson.gui.main_window import run_app

    run_app()


if __name__ == "__main__":
    main()
