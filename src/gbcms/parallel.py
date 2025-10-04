"""Parallel processing with joblib and Ray support."""

import logging
from typing import List, Callable, Any, Optional
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import os

from joblib import Parallel, delayed
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

logger = logging.getLogger(__name__)

# Try to import Ray (optional dependency)
try:
    import ray

    RAY_AVAILABLE = True
except ImportError:
    RAY_AVAILABLE = False
    logger.debug("Ray not available, using joblib for parallelization")


class ParallelProcessor:
    """Unified interface for parallel processing with joblib or Ray."""

    def __init__(
        self,
        n_jobs: int = -1,
        backend: str = "joblib",
        use_ray: bool = False,
        verbose: int = 0,
    ):
        """
        Initialize parallel processor.

        Args:
            n_jobs: Number of parallel jobs (-1 for all CPUs)
            backend: Backend to use ('joblib', 'ray', 'threading', 'multiprocessing')
            use_ray: Force use of Ray if available
            verbose: Verbosity level
        """
        self.n_jobs = n_jobs if n_jobs > 0 else os.cpu_count()
        self.backend = backend
        self.use_ray = use_ray and RAY_AVAILABLE
        self.verbose = verbose

        if self.use_ray and not RAY_AVAILABLE:
            logger.warning("Ray requested but not available, falling back to joblib")
            self.use_ray = False

        # Initialize Ray if requested
        if self.use_ray and not ray.is_initialized():
            ray.init(num_cpus=self.n_jobs, ignore_reinit_error=True)
            logger.info(f"Initialized Ray with {self.n_jobs} CPUs")

    def map(
        self,
        func: Callable,
        items: List[Any],
        description: str = "Processing",
        show_progress: bool = True,
    ) -> List[Any]:
        """
        Map function over items in parallel.

        Args:
            func: Function to apply
            items: Items to process
            description: Description for progress bar
            show_progress: Whether to show progress bar

        Returns:
            List of results
        """
        if self.use_ray:
            return self._map_ray(func, items, description, show_progress)
        else:
            return self._map_joblib(func, items, description, show_progress)

    def _map_joblib(
        self,
        func: Callable,
        items: List[Any],
        description: str,
        show_progress: bool,
    ) -> List[Any]:
        """Map using joblib."""
        if show_progress:
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TaskProgressColumn(),
            ) as progress:
                task = progress.add_task(f"[cyan]{description}...", total=len(items))

                results = []
                with Parallel(n_jobs=self.n_jobs, backend=self.backend) as parallel:
                    for result in parallel(delayed(func)(item) for item in items):
                        results.append(result)
                        progress.update(task, advance=1)

                return results
        else:
            return Parallel(n_jobs=self.n_jobs, backend=self.backend)(
                delayed(func)(item) for item in items
            )

    def _map_ray(
        self,
        func: Callable,
        items: List[Any],
        description: str,
        show_progress: bool,
    ) -> List[Any]:
        """Map using Ray."""
        # Convert function to Ray remote function
        remote_func = ray.remote(func)

        if show_progress:
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TaskProgressColumn(),
            ) as progress:
                task = progress.add_task(f"[cyan]{description}...", total=len(items))

                # Submit all tasks
                futures = [remote_func.remote(item) for item in items]

                # Collect results with progress
                results = []
                while futures:
                    # Wait for at least one task to complete
                    done, futures = ray.wait(futures, num_returns=1)
                    results.extend(ray.get(done))
                    progress.update(task, advance=len(done))

                return results
        else:
            futures = [remote_func.remote(item) for item in items]
            return ray.get(futures)

    def starmap(
        self,
        func: Callable,
        items: List[tuple],
        description: str = "Processing",
        show_progress: bool = True,
    ) -> List[Any]:
        """
        Starmap function over items (unpack tuples as arguments).

        Args:
            func: Function to apply
            items: List of tuples to unpack as arguments
            description: Description for progress bar
            show_progress: Whether to show progress bar

        Returns:
            List of results
        """
        if self.use_ray:
            return self._starmap_ray(func, items, description, show_progress)
        else:
            return self._starmap_joblib(func, items, description, show_progress)

    def _starmap_joblib(
        self,
        func: Callable,
        items: List[tuple],
        description: str,
        show_progress: bool,
    ) -> List[Any]:
        """Starmap using joblib."""
        wrapper = lambda args: func(*args)
        return self._map_joblib(wrapper, items, description, show_progress)

    def _starmap_ray(
        self,
        func: Callable,
        items: List[tuple],
        description: str,
        show_progress: bool,
    ) -> List[Any]:
        """Starmap using Ray."""
        wrapper = lambda args: func(*args)
        return self._map_ray(wrapper, items, description, show_progress)

    def shutdown(self):
        """Shutdown parallel processing resources."""
        if self.use_ray and ray.is_initialized():
            ray.shutdown()
            logger.info("Ray shutdown complete")


class BatchProcessor:
    """Process items in batches for better performance."""

    def __init__(
        self,
        batch_size: int = 100,
        n_jobs: int = -1,
        backend: str = "joblib",
    ):
        """
        Initialize batch processor.

        Args:
            batch_size: Number of items per batch
            n_jobs: Number of parallel jobs
            backend: Parallelization backend
        """
        self.batch_size = batch_size
        self.processor = ParallelProcessor(n_jobs=n_jobs, backend=backend)

    def process_batches(
        self,
        func: Callable,
        items: List[Any],
        description: str = "Processing batches",
    ) -> List[Any]:
        """
        Process items in batches.

        Args:
            func: Function to apply to each batch
            items: Items to process
            description: Description for progress

        Returns:
            Flattened list of results
        """
        # Create batches
        batches = [items[i : i + self.batch_size] for i in range(0, len(items), self.batch_size)]

        logger.info(f"Processing {len(items)} items in {len(batches)} batches")

        # Process batches in parallel
        batch_results = self.processor.map(func, batches, description)

        # Flatten results
        results = []
        for batch_result in batch_results:
            if isinstance(batch_result, list):
                results.extend(batch_result)
            else:
                results.append(batch_result)

        return results

    def shutdown(self):
        """Shutdown processor."""
        self.processor.shutdown()


def parallel_map(
    func: Callable,
    items: List[Any],
    n_jobs: int = -1,
    backend: str = "joblib",
    description: str = "Processing",
    show_progress: bool = True,
) -> List[Any]:
    """
    Convenience function for parallel mapping.

    Args:
        func: Function to apply
        items: Items to process
        n_jobs: Number of parallel jobs
        backend: Backend to use
        description: Progress description
        show_progress: Show progress bar

    Returns:
        List of results
    """
    processor = ParallelProcessor(n_jobs=n_jobs, backend=backend)
    try:
        return processor.map(func, items, description, show_progress)
    finally:
        processor.shutdown()


def parallel_starmap(
    func: Callable,
    items: List[tuple],
    n_jobs: int = -1,
    backend: str = "joblib",
    description: str = "Processing",
    show_progress: bool = True,
) -> List[Any]:
    """
    Convenience function for parallel starmapping.

    Args:
        func: Function to apply
        items: List of argument tuples
        n_jobs: Number of parallel jobs
        backend: Backend to use
        description: Progress description
        show_progress: Show progress bar

    Returns:
        List of results
    """
    processor = ParallelProcessor(n_jobs=n_jobs, backend=backend)
    try:
        return processor.starmap(func, items, description, show_progress)
    finally:
        processor.shutdown()


# Ray-specific utilities
if RAY_AVAILABLE:

    @ray.remote
    class VariantCounterActor:
        """Ray actor for stateful variant counting."""

        def __init__(self, config_dict: dict):
            """Initialize counter with configuration."""
            from .models import GetBaseCountsConfig

            self.config = GetBaseCountsConfig(**config_dict)
            self.processed_count = 0

        def count_variant_block(self, block_data: dict) -> dict:
            """
            Count a block of variants.

            Args:
                block_data: Dictionary with variant block information

            Returns:
                Dictionary with counting results
            """
            # Import here to avoid circular dependencies
            from .counter import BaseCounter

            counter = BaseCounter(self.config)
            # Process block
            # ... counting logic ...

            self.processed_count += 1
            return {"status": "success", "processed": self.processed_count}

        def get_stats(self) -> dict:
            """Get processing statistics."""
            return {"processed_blocks": self.processed_count}

    def create_ray_actors(n_actors: int, config_dict: dict) -> List:
        """
        Create Ray actors for distributed processing.

        Args:
            n_actors: Number of actors to create
            config_dict: Configuration dictionary

        Returns:
            List of actor handles
        """
        if not ray.is_initialized():
            ray.init(ignore_reinit_error=True)

        actors = [VariantCounterActor.remote(config_dict) for _ in range(n_actors)]

        logger.info(f"Created {n_actors} Ray actors")
        return actors

    def distribute_work_to_actors(
        actors: List,
        work_items: List[Any],
        description: str = "Distributing work",
    ) -> List[Any]:
        """
        Distribute work items to Ray actors.

        Args:
            actors: List of Ray actor handles
            work_items: Work items to process
            description: Progress description

        Returns:
            List of results
        """
        n_actors = len(actors)

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
        ) as progress:
            task = progress.add_task(f"[cyan]{description}...", total=len(work_items))

            # Distribute work round-robin
            futures = []
            for i, item in enumerate(work_items):
                actor = actors[i % n_actors]
                future = actor.count_variant_block.remote(item)
                futures.append(future)

            # Collect results
            results = []
            while futures:
                done, futures = ray.wait(futures, num_returns=1)
                results.extend(ray.get(done))
                progress.update(task, advance=len(done))

        return results
