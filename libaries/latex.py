
from libaries import general as gl
import os
import subprocess
import time
from threading import Thread
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


def compile_lualatex(tex_file, pdf_path=None, miktex_lualatex_path='C:/temp/MikTex/miktex/bin/x64/lualatex.exe', biber_path='C:/temp/MikTex/miktex/bin/x64/biber.exe'):
    def press_return_periodically(process):
        # Function to send newline every second to simulate "return"
        while process.poll() is None:  # Run as long as the process is active
            process.stdin.write('\n')
            process.stdin.flush()
            time.sleep(1)

    def run_subprocess(command, cwd):
        with subprocess.Popen(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, stdin=subprocess.PIPE, text=True, cwd=cwd) as process:
            # Start a separate thread to press "return" every second
            thread = Thread(target=press_return_periodically, args=(process,))
            thread.start()
            process.wait()  # Wait for the process to complete
            thread.join()  # Ensure the thread has finished

    run_path = os.path.dirname(os.path.realpath(__file__))
    txt_path = os.path.dirname(tex_file)
    main_name = os.path.basename(tex_file).removesuffix('.tex')

    output_pdf = pdf_path if pdf_path else os.path.join(txt_path, f"{main_name}.pdf")

    # Compile the .tex file using LuaLaTeX
    try:
        run_subprocess([miktex_lualatex_path, tex_file, '-output-directory', txt_path], cwd=txt_path)
    except Exception as e:
        print(f"LuaLaTeX compilation failed: {e}")

    # Run the bibliography tool (Biber)
    for i in range(1):
        try:
            run_subprocess([biber_path, main_name, '-output-directory', txt_path], cwd=txt_path)
        except Exception as e:
            print(f"Biber execution failed: {e}")

    # Run LuaLaTeX two more times for references update
    for i in range(3):
        try:
            run_subprocess([miktex_lualatex_path, tex_file, '-output-directory', txt_path], cwd=txt_path)
        except Exception as e:
            print(f"LuaLaTeX pass {i + 2} failed: {e}")

    # Check if the output PDF was created
    if os.path.exists(output_pdf):
        print(f"PDF successfully created at: {output_pdf}")
        return output_pdf
    else:
        print("PDF was not created, something went wrong.")
        return None


def include_Fig(string, FigInfo):
    figure_template = ("\\begin{figure}[H] \n "
                       "\\includegraphics[width=?FIGURE_WIDTH\\textwidth]{?FIGURE_PATH} \n "
                       "\\caption{ \\textit{?CAPTION}} \n "
                       "\\label{fig:?FIGURE_NAME} \n"
                       "\\end{figure}")

    if FigInfo is not None:
        figure_latex = gl.alias(figure_template,
                                {"1": "?FIGURE_PATH",
                                 "2": "?CAPTION",
                                 "3": "?FIGURE_NAME",
                                 "4": "?FIGURE_WIDTH"},
                                {"1": FigInfo["path"],
                                 "2": FigInfo["caption"],
                                 "3": FigInfo.name,
                                 "4": f"{FigInfo['width']}"})

    else:
        figure_latex = '\n \\ \\ \\ no data available  \\ \\ \\ \n'

    lines = find_keyword(string, '?FIG')
    string_out, _ = include_str(string, figure_latex, line=lines[0], replace=True)

    return string_out


def include_MultiFig(string, FigInfo):

    figure_template = ("\\begin{figure}[H] \n "
                       "\\includegraphics[width=?FIGURE_WIDTH\\textwidth]{?FIGURE_PATH} \n "
                       "\\caption{ \\textit{?CAPTION}} \n "
                       "\\label{fig:?FIGURE_NAME} \n"
                       "\\end{figure}")

    temp = []

    FigInfo = [info for info in FigInfo if info is not None]
    fig_string = ""

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
                       "\\captionsetup{type=table} \n"
                       "\\caption{ \\textit{?CAPTION}} \n "
                       "\\includegraphics[width=?FIGURE_WIDTH\\textwidth ]{?FIGURE_PATH} \n "
                       "\\label{tab:?FIGURE_NAME} \n"
                       "\\end{figure}")

    figure_template_no_cap = ("\\begin{figure}[H] \n "
                       "\\includegraphics[width=?FIGURE_WIDTH\\textwidth ]{?FIGURE_PATH} \n "
                       "\\label{tab:?FIGURE_NAME} \n"
                       "\\end{figure}")

    if FigInfo is not None:

        if FigInfo["caption"] is None:

            figure_latex = gl.alias(figure_template_no_cap,
                                    {"1": "?FIGURE_PATH",
                                     "3": "?FIGURE_NAME",
                                     "4": "?FIGURE_WIDTH"},
                                    {"1": FigInfo["path"],
                                     "3": FigInfo.name,
                                     "4": f"{FigInfo['width']}"})
        else:

            figure_latex = gl.alias(figure_template,
                                    {"1": "?FIGURE_PATH",
                                     "2": "?CAPTION",
                                     "3": "?FIGURE_NAME",
                                     "4": "?FIGURE_WIDTH"},
                                    {"1": FigInfo["path"],
                                     "2": FigInfo["caption"],
                                     "3": FigInfo.name,
                                     "4": f"{FigInfo['width']}"})

    else:
        figure_latex = '\\'

    lines = find_keyword(string, '?TABLE')
    string_out, _ = include_str(string, figure_latex, line=lines[0], replace=True)

    return string_out
