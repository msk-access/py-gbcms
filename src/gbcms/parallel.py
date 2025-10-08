"""Parallel processing with joblib backend only."""

import logging
import os
from collections.abc import Callable
from typing import Any

from joblib import Parallel, delayed
from rich.progress import BarColumn, Progress, TaskProgressColumn, TextColumn

logger = logging.getLogger(__name__)


class ParallelProcessor:
    """Unified interface for parallel processing with joblib."""

    def __init__(
        self,
        n_jobs: int = -1,
        backend: str = "joblib",
        verbose: int = 0,
    ):
        """
        Initialize parallel processor.

        Args:
            n_jobs: Number of parallel jobs (-1 for all CPUs)
            backend: Backend to use ('joblib', 'threading', 'multiprocessing')
            verbose: Verbosity level
        """
        self.n_jobs = n_jobs if n_jobs > 0 else os.cpu_count()
        self.backend = backend
        self.verbose = verbose

        # Map user-friendly backend names to joblib backends
        backend_map = {
            "joblib": "loky",  # Robust joblib backend
            "threading": "threading",  # Pure threading
            "multiprocessing": "multiprocessing",  # Process-based
            "loky": "loky",  # Explicit loky
        }
        self.joblib_backend = backend_map.get(backend, "loky")

        logger.debug(
            f"Initialized parallel processor with {self.n_jobs} jobs using {backend} -> {self.joblib_backend} backend"
        )

    def map(
        self,
        func: Callable,
        items: list,
        description: str = "Processing",
        show_progress: bool = True,
    ) -> list[Any]:
        """
        Apply function to each item in parallel.

        Args:
            func: Function to apply
            items: List of items to process
            description: Progress description
            show_progress: Whether to show progress bar

        Returns:
            List of results
        """
        return self._map_joblib(func, items, description, show_progress)

    def _map_joblib(
        self,
        func: Callable,
        items: list,
        description: str = "Processing",
        show_progress: bool = True,
    ) -> list[Any]:
        """Map using joblib."""
        if show_progress and len(items) > 10:
            # Use progress bar for larger workloads
            progress_columns = [
                TextColumn("[bold blue]{task.description}"),
                BarColumn(),
                TaskProgressColumn(),
                TextColumn("({task.completed}/{task.total})"),
            ]

            with Progress(*progress_columns, refresh_per_second=10) as progress:
                task = progress.add_task(description, total=len(items))

                def progress_wrapper(item):
                    result = func(item)
                    progress.update(task, advance=1)
                    return result

                with Parallel(n_jobs=self.n_jobs, backend=self.joblib_backend) as parallel:
                    return list(parallel(delayed(progress_wrapper)(item) for item in items))
        else:
            # Simple parallel execution
            with Parallel(n_jobs=self.n_jobs, backend=self.joblib_backend) as parallel:
                return list(parallel(delayed(func)(item) for item in items))

    def starmap(
        self,
        func: Callable,
        items: list,
        description: str = "Processing",
        show_progress: bool = True,
    ) -> list[Any]:
        """
        Apply function with arguments to each item in parallel.

        Args:
            func: Function to apply
            items: List of argument tuples
            description: Progress description
            show_progress: Whether to show progress bar

        Returns:
            List of results
        """

        def wrapper(args):
            return func(*args)

        return self._map_joblib(wrapper, items, description, show_progress)

    def shutdown(self):
        """Shutdown parallel processing resources."""
        # joblib handles cleanup automatically
        pass
