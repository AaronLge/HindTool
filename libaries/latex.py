import pandas as pd
import os
import subprocess
import time
from threading import Thread

from libaries import general as gl
from libaries import hindtoolplot as hc_plt
import shutil

def insertLatexVars(string, replacements):
    """
    Wrapper for the 'gl.alias' function to replace variables in the input string with specified replacements.

    Args:
        string (str): The input string containing variables to be replaced.
        replacements (dict): A dictionary where keys are the variables to be replaced
                              and values are their corresponding replacement values.

    Returns:
        str: A new string with the variables replaced by the provided replacements.

    Example:
        input_string = "This is a test with ?var1 and ?var2."
        replacements = {"var1": "value1", "var2": "value2"}
        output_string = insertLatexVars(input_string, replacements)
        # output_string will be "This is a test with value1 and value2."
    """

    Dummy_Orig = {Replacement_col: "?" + Replacement_col for Replacement_col in replacements.keys()}

    out_string = gl.alias(string, Dummy_Orig, replacements)

    return out_string


def find_keyword(string, keyword):
    """
    Finds the indices of lines in the input string that contain the specified keyword.

    Args:
        string (str): The input string to search within.
        keyword (str): The keyword to search for.

    Returns:
        list: A list of indices of the lines that contain the keyword.

    Example:
        input_string = "Line 1\nLine 2 with keyword\nLine 3"
        keyword = "keyword"
        indices = find_keyword(input_string, keyword)
        # indices will be [1]
    """
    indizes = []
    lines = string.split('\n')
    for i, line in enumerate(lines):
        if keyword in line:
            indizes.append(i)

    return indizes


def include_include(main, Include, line=None):
    """
    Inserts an `\\include{}` statement for a specified LaTeX file into the main LaTeX content.

    Parameters:
    -----------
    main : str
        The main LaTeX content as a string.
    Include : str
        The name of the LaTeX file (without the `.tex` extension) to include.
    line : int, optional
        The line number where the `\\include{}` statement should be inserted.
        If None, the statement is inserted one line above the `\\end{document}` line.

    Returns:
    --------
    tuple
        - `updated_main` (str): The modified LaTeX content with the inclusion statement added.
        - `insert_idx` (int): The index where the inclusion statement was inserted.
    """

    lines = main.split('\n')
    insert_idx = None
    if line is None:
        # Find the positions of "\begin{document}", "\end{document}", and the last "\include{"
        for i, line in enumerate(lines):
            if "\\end{document}" in line:
                end_doc_idx = i
                insert_idx = end_doc_idx - 1

        # Check "\end{document}" exist
        if insert_idx is None:
            raise ValueError("The file does not contain a valid LaTeX document structure.")

    else:
        insert_idx = line
    # Insert the string
    lines.insert(insert_idx, '\\include{' + f"{Include}" + "}")

    # Write the modified contents back to the file
    return '\n'.join(lines), insert_idx


def include_str(main, string, line, replace=False):
    """
    Inserts or replaces a string in the specified line of the main LaTeX content.

    Parameters:
    -----------
    main : str
        The main LaTeX content as a string.
    string : str
        The string to be inserted into the LaTeX content.
    line : int
        The line number where the string should be inserted.
    replace : bool, optional
        If True, replaces the existing content at the specified line.
        If False, the string is inserted without replacing existing content. Default is False.
    """


    lines_include = len(string.split('\n'))

    # Read the file content
    lines = main.split('\n')

    # Insert the string
    if replace:
        lines.pop(line)

    lines.insert(line, string)

    return '\n'.join(lines), line + lines_include


def compile_lualatex(tex_file, pdf_path=None, miktex_lualatex_path='C:/temp/MikTex/miktex/bin/x64/lualatex.exe', biber_path='C:/temp/MikTex/miktex/bin/x64/biber.exe'):
    """
    Compiles a LaTeX document using LuaLaTeX and Biber, ensuring all references and bibliographies are processed.

    Parameters:
    -----------
    tex_file : str
        The path to the `.tex` file to compile.
    pdf_path : str, optional
        The desired path for the output PDF. If not provided, the PDF will be saved in the same directory as `tex_file`.
    miktex_lualatex_path : str, optional
        Path to the LuaLaTeX executable. Defaults to 'C:/temp/MikTex/miktex/bin/x64/lualatex.exe'.
    biber_path : str, optional
        Path to the Biber executable. Defaults to 'C:/temp/MikTex/miktex/bin/x64/biber.exe'.

    Returns:
    --------
    str or None
        The path to the generated PDF if successful, or `None` if the PDF generation failed.

    Notes:
    ------
    - LuaLaTeX is run multiple times: initially to process the document, and two additional times to update references.
    - Biber is used to handle the bibliography, and it is executed between LuaLaTeX passes.
    - A separate thread sends newline characters to the LuaLaTeX process periodically, which may be necessary to handle specific interactive prompts.

    """
    def press_return_periodically(process):
        """
        Sends a newline character (`\\n`) to the input of a subprocess periodically.

        Parameters:
        -----------
        process : subprocess.Popen
            The process to which the newline character will be sent.
        """
        while process.poll() is None:  # Check if the process is active
            process.stdin.write('\n')
            process.stdin.flush()
            time.sleep(1)

    def run_subprocess(command, cwd):
        """
        Runs a subprocess command and handles periodic interaction with the process.

        Parameters:
        -----------
        command : list of str
            The command to execute as a list of arguments.
        cwd : str
            The working directory where the command should be executed.
        """
        with subprocess.Popen(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, stdin=subprocess.PIPE, text=True, cwd=cwd) as process:
            thread = Thread(target=press_return_periodically, args=(process,))
            thread.start()
            process.wait()  # Wait for the process to complete
            thread.join()  # Ensure the thread finishes

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

    # Run LuaLaTeX two more times to update references
    for i in range(3):
        try:
            run_subprocess([miktex_lualatex_path, tex_file, '-output-directory', txt_path], cwd=txt_path)
        except Exception as e:
            print(f"LuaLaTeX pass {i + 2} failed: {e}")

    # Verify if the output PDF was created
    if os.path.exists(output_pdf):
        print(f"PDF successfully created at: {output_pdf}")
        return output_pdf
    else:
        print("PDF was not created, something went wrong.")
        return None


def include_Fig(string, FigInfo):
    """
        Inserts a LaTeX figure into a string at a specific placeholder.

        Parameters:
        -----------
        string : str
            The LaTeX string where the figure should be inserted.
        FigInfo : pandas.Series or None
            Information about the figure to be included. Must contain the keys:
            - "path": The file path to the figure image.
            - "caption": The caption for the figure.
            - "width": The width of the figure (relative to `\textwidth`).
            - "name": The label for the figure (used for references).

            If `FigInfo` is `None`, a placeholder message is inserted.

        Returns:
        --------
        str
            The updated LaTeX string with the figure included.

        Notes:
        ------
        - The placeholder `?FIG` in the `string` marks where the figure will be inserted.
        - If `FigInfo` is not provided, a default message (`No data available.`) will replace the placeholder.
        """

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
        figure_latex = 'No data available. \n'

    lines = find_keyword(string, '?FIG')
    string_out, _ = include_str(string, figure_latex, line=lines[0], replace=True)

    return string_out


def include_MultiFig(string, FigInfo):
    """
    Inserts multiple LaTeX figures into a string at a specific placeholder.

    Parameters:
    -----------
    string : str
        The LaTeX string where the figures should be inserted.
    FigInfo : list of pandas.Series or None
        A list containing information about each figure to be included. Each element must contain the keys:
        - "path": The file path to the figure image.
        - "caption": The caption for the figure.
        - "width": The width of the figure (relative to `\textwidth`).
        - "name": The label for the figure (used for references).

        If an element in the list is `None`, it is skipped.

    Returns:
    --------
    str
        The updated LaTeX string with all figures included.

    Notes:
    ------
    - The placeholder `?MULTIFIG` in the `string` marks where the figures will be inserted.
    - Figures are formatted using a common LaTeX `figure` environment template.
    """

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


def include_MultiTab(string, FigInfo):
    """
    Inserts multiple LaTeX table-style figures into a string at a specific placeholder.

    Parameters:
    -----------
    string : str
        The LaTeX string where the table-style figures should be inserted.
    FigInfo : list of pandas.Series or None
        A list containing information about each figure to be included. Each element must contain the keys:
        - "path": The file path to the figure image.
        - "caption": The caption for the figure.
        - "width": The width of the figure (relative to `\textwidth`).
        - "name": The label for the figure (used for references).

        If an element in the list is `None`, it is skipped.

    Returns:
    --------
    str
        The updated LaTeX string with all table-style figures included.

    Notes:
    ------
    - The placeholder `?MULTITAB` in the `string` marks where the figures will be inserted.
    - Figures are formatted to appear as tables using `\captionsetup{type=table}` in the LaTeX template.
     """

    figure_template = ("\\begin{figure}[H] \n "
                       "\\captionsetup{type=table} \n"
                       "\\caption{ \\textit{?CAPTION}} \n "
                       "\\includegraphics[width=?FIGURE_WIDTH\\textwidth ]{?FIGURE_PATH} \n "
                       "\\label{tab:?FIGURE_NAME} \n"
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

    lines = find_keyword(string, '?MULTITAB')
    string_out, _ = include_str(string, fig_string, line=lines[0], replace=True)

    return string_out


def include_TableFig(string, FigInfo):
    """
    Inserts a single table-style figure into a LaTeX string at a specific placeholder.

    Parameters:
    -----------
    string : str
        The LaTeX string where the table-style figure should be inserted.
    FigInfo : dict or None
        A dictionary containing information about the figure to be included. Must include the keys:
        - "path": The file path to the figure image.
        - "caption": The caption for the figure (or None for no caption).
        - "width": The width of the figure (relative to `\textwidth`).
        - "name": The label for the figure (used for references).

        If `FigInfo` is `None`, the function inserts an empty placeholder (`\\`).

    Returns:
    --------
    str
        The updated LaTeX string with the table-style figure included.

    Notes:
    ------
    - The placeholder `?TABLE` in the `string` marks where the figure will be inserted.
    - Figures are formatted to appear as tables using `\captionsetup{type=table}` in the LaTeX template.
    """

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


def initilize_document(DocumentMeta, Revisions, bib_paths, acronyms_path, save_path, map=None, introduction_text=None, document_purpose_text=None):
    """
    Initializes a JBO LaTeX document by preparing metadata, bibliographies, acronyms, and other content for rendering.

    Parameters:
    -----------
    DocumentMeta : dict
        Metadata for the document, including keys like "RevisionJBO", "RevisionEmployer", and "RevisionDate".
        Values set to 'auto' will be replaced with corresponding values from the `Revisions` DataFrame.
    Revisions : pandas.DataFrame
        A table of revision data, where the last row provides default values for metadata fields.
    bib_paths : list of str
        Paths to the bibliography files to be included in the document.
    acronyms_path : str
        Path to the file containing acronyms, to be inserted into the document.
    save_path : str
        Directory where generated figures and files will be saved.
    map : tuple of str, optional
        A tuple containing the path and caption of a map image to include in the document. Default is None.
    introduction_text : str, optional
        Text to be included in the introduction section. Default is None.
    document_purpose_text : str, optional
        Text describing the purpose of the document, included in the introduction. Default is None.

    Returns:
    --------
    tuple
        - `main_tex` (str): The main LaTeX content with placeholders replaced.
        - `titlepage_tex` (str): LaTeX content for the title page.
        - `introduction_tex` (str): LaTeX content for the introduction, including figures and tables.

    Raises:
    -------
    Any exception raised by file operations or external functions will propagate to the caller.
    """

    curr_path = os.path.dirname(os.path.realpath(__file__))
    figsize_fullpage = [17 * 0.39370079, None]

    if DocumentMeta["RevisionJBO"] == 'auto':
        DocumentMeta["RevisionJBO"] = Revisions.iloc[-1, 0]

    if DocumentMeta["RevisionEmployer"] == 'auto':
        DocumentMeta["RevisionEmployer"] = Revisions.iloc[-1, 1]

    if DocumentMeta["RevisionDate"] == 'auto':
        DocumentMeta["RevisionDate"] = Revisions.iloc[-1, 2]

    if DocumentMeta["RevisionDate"] == 'auto':
        DocumentMeta["RevisionDate"] = Revisions.iloc[-1, 2]

    main_path = curr_path + '\\latex_main_template.txt'
    titlepage_path = curr_path + '\\latex_titlepage_template.txt'
    introduction_path = curr_path + '\\latex_introduction_template.txt'

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

    shutil.copyfile(curr_path + "\\latex_Status_table.png", save_path + "\\Status_table.png")
    shutil.copyfile(curr_path+'\\JBO_logo.jpg', save_path + '\\JBO_logo.jpg')

    pic = "Status_table"
    FIGURES.loc[pic, "filename"] = f"{pic}.png"
    FIGURES.loc[pic, "path"] = save_path + f"\\{pic}.png"
    FIGURES.loc[pic, "caption"] = None
    FIGURES.loc[pic, "width"] = 0.4

    # Figures
    pic = "Revision_Table_page_1"
    FIGURES.loc[pic, "filename"] = f"{pic}.png"
    FIGURES.loc[pic, "path"] = save_path + f"\\{pic}.png"
    FIGURES.loc[pic, "caption"] = None
    FIGURES.loc[pic, "width"] = 1

    if map is not None:
        pic = "map"
        FIGURES.loc[pic, "filename"] = f"{pic}.png"
        FIGURES.loc[pic, "path"] = map[0]
        FIGURES.loc[pic, "caption"] = map[0]
        FIGURES.loc[pic, "width"] = 0.4

    FIGURES.loc[:, "path"] = [string.replace("\\", "/") for string in FIGURES.loc[:, "path"]]

    with open(main_path, 'r', encoding='utf-8') as file:
        main_tex = file.read()

    with open(titlepage_path, 'r', encoding='utf-8') as file:
        titlepage_tex = file.read()

    with open(introduction_path, 'r', encoding='utf-8') as file:
        introduction_tex = file.read()

    with open(acronyms_path, 'r', encoding='utf-8') as file:
        acronyms = file.read()

    main_tex = insertLatexVars(main_tex, {"ACRONYMS": acronyms})

    # import biblografys:
    biblografies = []
    for bib_path in bib_paths:
        bib_path = bib_path.replace('\\','/')
        biblografies.append("\\addbibresource{" + f"{bib_path}" + "}")

    biblografies = '\n'.join(biblografies)
    main_tex = insertLatexVars(main_tex, {'Biblografies': biblografies})

    # insert titlepage
    chapter = 'titlepage'
    titlepage_tex = insertLatexVars(titlepage_tex, DocumentMeta)
    main_tex, last_idx = include_include(main_tex, chapter)

    #pagestyle
    main_tex, last_idx = include_str(main_tex, '\\pagestyle{fancy}', last_idx + 1)
    main_tex = insertLatexVars(main_tex, DocumentMeta)

    # insert introduction
    chapter = 'introduction'

    introduction_tex = include_TableFig(introduction_tex, FIGURES.loc["Revision_Table_page_1"])
    introduction_tex = include_TableFig(introduction_tex, FIGURES.loc["Status_table"])
    introduction_tex = insertLatexVars(introduction_tex, {"IntroductionText": introduction_text, "DocumentPurposeText": document_purpose_text})

    introduction_tex = include_Fig(introduction_tex, FIGURES.loc["map"] if "map" in FIGURES.index else None)

    main_tex, last_idx = include_include(main_tex, chapter)

    return main_tex, titlepage_tex, introduction_tex
