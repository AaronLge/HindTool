import pandas as pd

from libaries import general as gl
from libaries import hindtoolplot as hc_plt

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
            if "\\end{document}" in line:
                end_doc_idx = i
                insert_idx = end_doc_idx - 1

        # Check "\end{document}" exist
        if end_doc_idx is None:
            raise ValueError("The file does not contain a valid LaTeX document structure.")

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
        figure_latex = '\n'

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


def initilize_document(DocumentMeta, Revisions, bib_paths, acronyms_path, save_path):

    figsize_fullpage = [17 * 0.39370079, None]


    if Revisions["DocumentStatus"] == 'auto':
        if Revisions.iloc[-1, 3] == 'First Release':
            Revisions["DocumentStatus"] = 'FR'
        elif Revisions.iloc[-1, 3] == 'Final':
            Revisions["DocumentStatus"] = 'FIN'
        elif Revisions.iloc[-1, 3] == 'Preliminary':
            Revisions["DocumentStatus"] = 'PRE'
        elif Revisions.iloc[-1, 3] == 'Issued for Review':
            Revisions["DocumentStatus"] = 'IFR'
        elif Revisions.iloc[-1, 3] == 'Issued for Excecution':
            Revisions["DocumentStatus"] = 'IFE'
        elif Revisions.iloc[-1, 3] == 'As Built':
            Revisions["DocumentStatus"] = 'AB'
        elif Revisions.iloc[-1, 3] == 'First Release':
            Revisions["DocumentStatus"] = 'AB'

        else:
            print(f"initilize document: {Revisions.iloc[-1, 3]} no known keyword, leaving it as is")
            Revisions["DocumentStatus"] = Revisions.iloc[-1, 3]

    if Revisions["RevisionJBO"] == 'auto':
        Revisions["RevisionJBO"] = Revisions.iloc[-1, 0]

    if Revisions["RevisionEmployer"] == 'auto':
        Revisions["RevisionEmployer"] = Revisions.iloc[-1, 1]

    if Revisions["RevisionDate"] == 'auto':
        Revisions["RevisionDate"] = Revisions.iloc[-1, 2]

    if Revisions["RevisionDate"] == 'auto':
        Revisions["RevisionDate"] = Revisions.iloc[-1, 2]

    main_path = 'latex_main_template.txt'
    titlepage_path = 'latex_titlepage_template.txt'
    introduction_path = 'latex_introduction_template.txt'

    FIGURES = pd.DataFrame()

    # Figures
    # Revision Table
    FIG = hc_plt.table(Revisions.values,
                       collabels=Revisions.columns,
                       rowlabels=None,
                       row_label_name='Parameters',
                       figsize=figsize_fullpage,
                       cell_height=0.7,
                       cell_width=[1, 1, 2, 3],
                       use_pgf=True)

    gl.save_figs_as_png([FIG], save_path + f'\\Revision_Table', dpi=500)

    pic = "latex_Status_table"
    FIGURES.loc[pic, "filename"] = f"{pic}.jpg"
    FIGURES.loc[pic, "path"] = save_path + f"\\{pic}.jpg"
    FIGURES.loc[pic, "caption"] = None
    FIGURES.loc[pic, "width"] = 0.4

    # Figures
    pic = "Revision_table"
    FIGURES.loc[pic, "filename"] = f"{pic}.jpg"
    FIGURES.loc[pic, "path"] = save_path + f"\\{pic}.jpg"
    FIGURES.loc[pic, "caption"] = None
    FIGURES.loc[pic, "width"] = 0.4

    with open(main_path, 'r', encoding='utf-8') as file:
        main_tex = file.read()

    with open(titlepage_path, 'r', encoding='utf-8') as file:
        titlepage_tex = file.read()

    with open(introduction_path, 'r', encoding='utf-8') as file:
        introduction_tex = file.read()

    with open(acronyms_path, 'r', encoding='utf-8') as file:
        acronyms_path = file.read()

    # inport acronyms
    with open(acronyms_path, 'r') as f:
        acros = f.read()
    main_tex = insertLatexVars(main_tex, {"ACRONYMS": acros})

    # import biblografys:
    biblografies = []
    for bib_path in bib_paths:
        bib_path = bib_path.replace('\\','/')
        biblografies.append("\\addbibresource{" + f"{bib_path}" + "}")

    biblografies = '\n'.join(biblografies)
    main_tex = insertLatexVars(main_tex, {'?Biblografies': biblografies})

    # insert titlepage
    chapter = 'titlepage'
    titlepage_tex = insertLatexVars(titlepage_tex, DocumentMeta)
    main_tex, last_idx = include_include(main_tex, chapter)

    #pagestyle
    main_tex, last_idx = include_str(main_tex, '\\pagestyle{fancy}', last_idx + 1)

    # insert introductiot
    chapter = 'introduction'
    introduction_tex = insertLatexVars(introduction_tex, DocumentMeta)

    introduction_tex = include_TableFig(introduction_tex, FIGURES.loc["Revision_Table_page_1"])
    introduction_tex = include_TableFig(introduction_tex, FIGURES.loc["latex_Status_table"])

    main_tex, last_idx = include_include(main_tex, chapter)

    return main_tex, titlepage_tex, introduction_tex