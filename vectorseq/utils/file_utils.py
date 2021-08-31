from pathlib import Path
import pickle


def create_dir(filepath: Path):
    "create diretory if it doesn't exist in the filepath"
    # os.makedirs(os.path.dirname(filepath), exist_ok=True)
    Path(filepath).mkdir(parents=True, exist_ok=True)
    return filepath


def file_rename(save_directory: Path, original_name: str, new_name: str):
    """
    rename saved figure to something else
    """
    pname = save_directory / f"{original_name}"
    pname.rename(save_directory / f"{new_name}")


def pickling(data_to_pickle, filename: str):
    """
    pickle save files
    """
    with open(filename, "wb+") as f:
        pickle.dump(data_to_pickle, f)


def unpickling(filename: str):
    """
    for unpickling files
    """
    with open(filename, "rb+") as f:
        unpickled_file = pickle.load(f)
    return unpickled_file
