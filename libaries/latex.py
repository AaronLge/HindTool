import sys
path = r"C:\\temp\\python_self_crated\\packages"
sys.path.insert(0, path)
from allib import general as gl


def insertLatexVars(path_template, Replacements):

    with open(path_template, 'r', encoding='utf-8') as file:
        lines = file.read()

    Dummy_Orig = {Replacement_col: "?" + Replacement_col for Replacement_col in Replacements.keys()}

    out_string = gl.alias(lines, Dummy_Orig, Replacements)

    return out_string


def include_include(main, Include):
    # Read the file content

    # Initialize variables to track the position of key elements
    begin_doc_idx = None
    end_doc_idx = None
    last_include_idx = None

    lines = main.split('\n')

    # Find the positions of "\begin{document}", "\end{document}", and the last "\include{"
    for i, line in enumerate(lines):
        if "\\begin{document}" in line:
            begin_doc_idx = i
        if "\\end{document}" in line:
            end_doc_idx = i
        if line.strip().startswith("\\include{"):
            last_include_idx = i

    # Check if "\begin{document}" and "\end{document}" exist
    if begin_doc_idx is None or end_doc_idx is None:
        raise ValueError("The file does not contain a valid LaTeX document structure.")

    # Determine where to insert the string
    if last_include_idx is not None:
        # Insert after the last "\include{"
        insert_idx = last_include_idx + 1
    else:
        # Insert after "\begin{document}"
        insert_idx = begin_doc_idx + 1

    # Insert the string
    lines.insert(insert_idx, '\\include{' + f"{Include}" + "}")

    # Write the modified contents back to the file

    return '\n'.join(lines)