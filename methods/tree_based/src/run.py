import argparse
import sys

sys.path.append('config')
from model import test_model, train_model, tune_model

def main():
    parser = argparse.ArgumentParser(
        description="A command line tool to manage the project")
    parser.add_argument('stage',
                        metavar='stage',
                        type=str,
                        choices=['tune', 'train', 'test'],
                        help="Stage to run.")

    stage = parser.parse_args().stage

    if stage == "tune":
        print("Tuning model...")
        tune_model()

    if stage == "train":
        print("Train model...")
        train_model()

    elif stage == "test":
        print("Test model...")
        test_model()

if __name__ == "__main__":
    main()
