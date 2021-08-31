from tqdm.auto import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from typing import Iterable


def parallel_process(
    iterable: Iterable,
    function: callable,
    n_jobs: int = cpu_count(),
    use_args: bool = False,
    use_kwargs: bool = False,
    desc: str = "",
):
    """
    A parallel version of the map function with a progress bar.

    Args:
        iterable (collection): An array-like or dict-like to iterate over
        function (function): A python function to apply to the elements of array
        n_jobs (int, default=num cpu cores on machine): The number of cores to use
        use_args (boolean, default=False): Whether to consider the elements of array as tuples of arguments to function.
            Tuple elements are passed to function arguments by position.
        use_kwargs (boolean, default=False): Whether to consider the elements of array as dictionaries of
            keyword arguments to function
        desc (string, default=""): Description on progress bar
    Returns:
        [function(iterable[0]), function(iterable[1]), ...]
    """
    # If we set n_jobs to 1, just run a list comprehension. This is useful for benchmarking and debugging.
    if n_jobs == 1:
        if use_kwargs:
            return [function(**a) for a in tqdm(iterable)]
        elif use_args:
            return [function(*a) for a in tqdm(iterable)]
        else:
            return [function(a) for a in tqdm(iterable)]
    # Assemble the workers
    with ProcessPoolExecutor(max_workers=n_jobs) as pool:
        # Pass the elements of array into function
        if use_kwargs:
            futures = [pool.submit(function, **a) for a in iterable]
        elif use_args:
            futures = [pool.submit(function, *a) for a in iterable]
        else:
            futures = [pool.submit(function, a) for a in iterable]
        # Print out the progress as tasks complete
        kwargs = {
            "total": len(futures),
            "unit": "it",
            "unit_scale": True,
            "leave": True,
            "desc": f"{desc} (Dispatch)",
        }
        for f in tqdm(as_completed(futures), **kwargs):
            pass
    out = []
    # Get the results from the futures.
    for i, future in tqdm(enumerate(futures), desc=f"{desc} (Completed)"):
        try:
            out.append(future.result())
        except Exception as e:
            out.append(e, f"Occurred with input element at index: {i}.")
    return out
