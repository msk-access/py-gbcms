"""
Logging utilities for gbcms.

Provides centralized logging configuration with dual output:
- Structured logging via Python logging module
- Rich console output for interactive use
"""

import logging
import time
from collections.abc import Callable
from contextlib import contextmanager
from functools import wraps
from typing import Any

from rich.console import Console
from rich.logging import RichHandler

__all__ = [
    "setup_logging",
    "get_logger",
    "timed",
    "log_call",
]

# Module-level console for rich output
_console = Console()


def setup_logging(verbose: bool = False, log_file: str | None = None) -> None:
    """
    Configure logging for gbcms.

    Args:
        verbose: If True, set log level to DEBUG. Otherwise INFO.
        log_file: Optional path to write logs to file.
    """
    log_level = logging.DEBUG if verbose else logging.INFO

    handlers: list[logging.Handler] = [
        RichHandler(
            console=_console,
            rich_tracebacks=True,
            markup=True,
            show_path=verbose,
        )
    ]

    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        )
        handlers.append(file_handler)

    logging.basicConfig(
        level=log_level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=handlers,
        force=True,  # Override existing config
    )


def get_logger(name: str) -> logging.Logger:
    """Get a logger instance for the given module name."""
    return logging.getLogger(name)


@contextmanager
def timed(operation: str, logger: logging.Logger | None = None):
    """
    Context manager for timing operations.

    Args:
        operation: Description of the operation being timed.
        logger: Logger to use. If None, uses root logger.

    Example:
        with timed("Loading variants", logger):
            variants = load_variants()
    """
    log = logger or logging.getLogger(__name__)
    start = time.perf_counter()
    log.debug("Starting: %s", operation)
    try:
        yield
    finally:
        elapsed = time.perf_counter() - start
        log.debug("Completed: %s (%.3fs)", operation, elapsed)


def log_call(logger: logging.Logger | None = None) -> Callable:
    """
    Decorator to log function calls with timing.

    Args:
        logger: Logger to use. If None, uses function's module logger.

    Example:
        @log_call()
        def process_sample(sample_name: str) -> dict:
            ...
    """

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            log = logger or logging.getLogger(func.__module__)
            log.debug("Calling %s", func.__name__)
            start = time.perf_counter()
            try:
                result = func(*args, **kwargs)
                elapsed = time.perf_counter() - start
                log.debug("%s completed (%.3fs)", func.__name__, elapsed)
                return result
            except Exception as e:
                log.error("%s failed: %s", func.__name__, e)
                raise

        return wrapper

    return decorator
