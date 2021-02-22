"""Generated latex code from a CSV file. Can also escape special characters
that latex cannot support in a string by themselves, such as _."""
from pathlib import Path


def main(csv_path, sep=",", chars_to_escape=["_"]):
    csv_path = Path(csv_path)
    latex = "\\begin{table}\n\t\\centering\n\t\\begin{tabular}{|"
    with csv_path.open() as csv:
        header = next(csv).split(sep)
        for i in range(0, len(header)):
            latex += "c|"
        latex += "}\n\t\t\\hline\n\t\t"
        for word in header:
            word = word.replace("\n", "")
            word = escape_chars(word, chars_to_escape)
            latex += f"\\textbf{{{word}}} & "
        latex = latex[:-2] + "\\\\\n\t\t\\hline\\hline\n\t\t"
        for line in csv:
            line = line.replace('"', "")
            for word in line.split(sep):
                word = word.replace("\n", "")
                word = escape_chars(word, chars_to_escape)
                latex += f"{word} & "
            latex = latex[:-2] + "\\\\\n\t\t\\hline\n\t\t"
        latex += "\t\\end{tabular}\n\t\\label{tab:}\n\t\\caption{}\n\\end{table}"

    out_path = csv_path.parent / (csv_path.stem + "_latexed.txt")
    with out_path.open("w+") as out:
        out.writelines(latex.replace('"', ""))
    return True


def escape_chars(string, characters):
    if type(characters) is str:
        characters = list(characters)
    for character in characters:
        string = string.replace(str(character), ("\\" + character))
    return string


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", help="Path to the csv file.", type=str)
    parser.add_argument("--sep", help="Line separator", default=",")
    parser.add_argument(
        "--escape", help="Chars to escape", nargs="*", default=["_", "&", "%"]
    )

    args = parser.parse_args()
    main(
        args.csv_path,
        args.sep,
        args.escape,
    )
