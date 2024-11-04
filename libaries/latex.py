
from libaries import general as gl
import subprocess
import os
import shutil
import pexpect


def insertLatexVars(string, replacements):
    """Wrapper for 'gl.alias' function to replace variables in string with replacements.#
        '?' + replacements.keys() as original
        replacements.values() as replacements"""

    Dummy_Orig = {Replacement_col: "?" + Replacement_col for Replacement_col in replacements.keys()}

    out_string = gl.alias(string, Dummy_Orig, replacements)

    return out_string


def find_keyword(string, keyword):
    """finds lines in string in which there is the keyword"""

    indizes = []
    lines = string.split('\n')
    for i, line in enumerate(lines):
        if keyword in line:
            indizes.append(i)

    return indizes


def include_include(main, Include, line=None):
    # Read the file content

    # Initialize variables to track the position of key elements
    begin_doc_idx = None
    end_doc_idx = None
    last_include_idx = None

    lines = main.split('\n')

    if line is None:
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

    else:
        insert_idx = line
    # Insert the string
    lines.insert(insert_idx, '\\include{' + f"{Include}" + "}")

    # Write the modified contents back to the file

    return '\n'.join(lines), insert_idx


def include_str(main, string, line, replace=False):
    lines_include = len(string.split('\n'))

    # Read the file content
    lines = main.split('\n')

    # Insert the string
    if replace:
        lines.pop(line)

    lines.insert(line, string)

    return '\n'.join(lines), line + lines_include


def compile_lualatex(tex_file, pdf_path=None, miktex_lualatex_path='C:/temp/MikTex/miktex/bin/x64/lualatex.exe', biber_path = 'C:/temp/MikTex/miktex/bin/x64/biber.exe'):
    # def run_subprocess():
    #     command = f"{miktex_lualatex_path} {tex_file} -output-directory {txt_path}"
    #     with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, text=True, cwd=txt_path) as process:
    #         try:
    #             # Read output line by line
    #             while True:
    #                 output = process.stdout.readline()
    #                 if output:
    #                     print(output.strip())  # Print the output
    #
    #                 # Check for the prompt and send an empty line
    #                 if '?' in output:
    #                     process.stdin.write('\n')
    #                     process.stdin.write('\n')
    #                     process.stdin.write('\n')
    #                     process.stdin.write('\n')
    #
    #                     process.stdin.flush()
    #
    #         except Exception as e:
    #             print(f"An error occurred: {e}")

    """
    Compile a LaTeX document using LuaLaTeX and save the output PDF to a specified location.

    Parameters:
    tex_files (str): Full path to the LaTeX (.tex) file you want to compile.
    output_pdf (str, optional): Full path to save the generated PDF. If not provided, the PDF will be saved
                                     in the same directory as the .tex file with the same base name.
    miktex_lualatex_path (str, optional): Full path to the LuaLaTeX executable (default is set to MiKTeX's standard location: 'C:/temp/MikTex/miktex/bin/x64/lualatex.exe').

    Returns:
    str: Path to the generated PDF if successful, or None if compilation fails.
    """
    run_path = os.path.dirname(os.path.realpath(__file__))
    txt_path = os.path.dirname(tex_file)

    Main_name = os.path.basename(tex_file)
    Main_name = Main_name.removesuffix('.tex')

    if pdf_path is None:
        output_pdf = txt_path + '\\' + Main_name + '.pdf'

    else:
        output_pdf = pdf_path
        shutil.move(run_path + '\\' + Main_name + '.pdf', output_pdf)
    # Step 2: Compile the .tex file using LuaLaTeX binary from MiKTeX
    try:
        # Extract the directory where the PDF should be output

        subprocess.run(
            [miktex_lualatex_path, tex_file, '-output-directory', txt_path],
            check=True,
            cwd=txt_path,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        #run_subprocess()

        # # Run the bibliography tool (e.g., Biber)
        # subprocess.run(
        #     [biblio_tool_path, Main_name],
        #     check=True,
        #     cwd=tex_dir,
        #     stdout=subprocess.PIPE,
        #     stderr=subprocess.PIPE,
        #     timeout=timeout
        #)
    except subprocess.CalledProcessError as e:
        print(f"LuaLaTeX compilation failed (maybe): {e}")
    try:
        subprocess.run(
            [biber_path, Main_name, '-output-directory', txt_path],
            check=True,
            cwd=txt_path,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        #run_subprocess()
    except subprocess.CalledProcessError as e:
        print(f"LuaLaTeX compilation failed (maybe): {e}")
    try:
        subprocess.run(
            [miktex_lualatex_path, tex_file, '-output-directory', txt_path],
            check=True,
            cwd=txt_path,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        #run_subprocess()

    except subprocess.CalledProcessError as e:
        print(f"LuaLaTeX compilation failed (maybe): {e}")

    # Step 3: Verify if the output PDF was created
    if os.path.exists(output_pdf):
        print(f"PDF successfully created at: {output_pdf}")
        return output_pdf
    else:
        print(f"PDF was not created, something went wrong.")
        return output_pdf


def include_Fig(string, FigInfo):
    figure_template = ("\\begin{figure}[H] \n "
                       "\\centering \n "
                       "\\includegraphics[width=?FIGURE_WIDTH\\textwidth]{?FIGURE_PATH} \n "
                       "\\caption{?CAPTION} \n "
                       "\\label{fig:?FIGURE_NAME} \n"
                       "\\end{figure}")

    figure_latex = gl.alias(figure_template,
                            {"1": "?FIGURE_PATH",
                             "2": "?CAPTION",
                             "3": "?FIGURE_NAME",
                             "4": "?FIGURE_WIDTH"},
                            {"1": FigInfo["path"],
                             "2": FigInfo["caption"],
                             "3": FigInfo.name,
                             "4": f"{FigInfo['width']}"})

    lines = find_keyword(string, '?FIG')
    string_out, _ = include_str(string, figure_latex, line=lines[0], replace=True)

    return string_out

def include_MultiFig(string, FigInfo):

    figure_template = ("\\begin{figure}[H] \n "
                       "\\centering \n "
                       "\\includegraphics[width=?FIGURE_WIDTH\\textwidth]{?FIGURE_PATH} \n "
                       "\\caption{?CAPTION} \n "
                       "\\label{fig:?FIGURE_NAME} \n"
                       "\\end{figure}")

    temp = []

    for Fig in FigInfo:
        figure_latex = gl.alias(figure_template,
                                {"1": "?FIGURE_PATH",
                                 "2": "?CAPTION",
                                 "3": "?FIGURE_NAME",
                                 "4": "?FIGURE_WIDTH"},
                                {"1": Fig["path"],
                                 "2": Fig["caption"],
                                 "3": Fig.name,
                                 "4": f"{Fig['width']}"})

        temp.append(figure_latex)
        fig_string = "\n".join(temp)

    lines = find_keyword(string, '?MULTIFIG')
    string_out, _ = include_str(string, fig_string, line=lines[0], replace=True)

    return string_out

def include_TableFig(string, FigInfo):
    figure_template = ("\\begin{figure}[H] \n "
                       "\\centering \n "
                       "\\includegraphics[width=?FIGURE_WIDTH\\textwidth ]{?FIGURE_PATH} \n "
                       "\\captionsetup{type=table} \n"
                       "\\caption{?CAPTION} \n "
                       "\\label{tab:?FIGURE_NAME} \n"
                       "\\end{figure}")

    figure_latex = gl.alias(figure_template,
                            {"1": "?FIGURE_PATH",
                             "2": "?CAPTION",
                             "3": "?FIGURE_NAME",
                             "4": "?FIGURE_WIDTH"},
                            {"1": FigInfo["path"],
                             "2": FigInfo["caption"],
                             "3": FigInfo.name,
                             "4": f"{FigInfo['width']}"})

    lines = find_keyword(string, '?TABLE')
    string_out, _ = include_str(string, figure_latex, line=lines[0], replace=True)

    return string_out
