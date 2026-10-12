import inspect
import os
import subprocess
import sys
import traceback
from pathlib import Path

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'GMAP'
copyright = '2026, KE van Adrichem'
author = 'KE van Adrichem'
release = '1.0'  # will be overwritten

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'numpydoc',
    'sphinx_design',
    'sphinx.ext.mathjax',
    'sphinx.ext.linkcode'
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Variables needed in functions -------------------------------------------

# the path of the base repository folder
repopath = Path(__file__).parent.parent.resolve()
sys.path.append(str(repopath))

# The name of the current version (in the built docs)
latest_version = os.environ.get("LATEST_VERSION", "stable")
source_branch = os.environ.get("SOURCE_BRANCH", "main")
if source_branch == "development":
    version_match = "dev"
    release = "dev"
else:
    version_match = ".".join(source_branch.split(".")[:2])
    release = source_branch

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_theme_options = {
    "logo": {
        "image_light": "Figures/GMAP_banner_light.png",
        "image_dark": "Figures/GMAP_banner_dark.png",
    },
    "header_links_before_dropdown": 4,
    "github_url": "https://github.com/lacourjansenlab/GMAP",
    "navbar_end": [
        "search-button",
        "theme-switcher",
        "version-switcher",
        "navbar-icon-links",
    ],
    "navbar_persistent": [],
    "switcher": {
        "json_url": "https://kimvana.github.io/GMAP_clone/stable/_static/switcher.json",
        "version_match": version_match,
    },
    "check_switcher": False,
    "show_version_warning_banner": True,
    "sticky_banners": True,
}
html_static_path = ['_static']
html_favicon = "Figures/GMAP_ico.ico"

suppress_warnings = ["docutils"]

html_js_files = [(
    None,
    {
        "body": """
        document.addEventListener("DOMContentLoaded", function() {
            const targetDiv = document.getElementById("central-dynamic-block");
            if (targetDiv) {
                // DOCUMENTATION_OPTIONS.URL_ROOT is an internal Sphinx variable 
                // that automatically resolves the correct relative path to your site's root.
                const rootPath = typeof DOCUMENTATION_OPTIONS !== 'undefined' ? DOCUMENTATION_OPTIONS.URL_ROOT : './';
                const targetUrl = `${rootPath}_static/central_block.html`;

                fetch(targetUrl, { cache: "no-store" })
                    .then(response => {
                        if (!response.ok) throw new Error("File not found");
                        return response.text();
                    })
                    .then(data => {
                        targetDiv.innerHTML = data;
                    })
                    .catch(error => {
                        console.error("Error updating live block:", error);
                        targetDiv.innerHTML = "<p style='color: red;'>Failed to load live updates.</p>";
                    });
            }}
        );
        """
    }
)]

# -- Other options -----------------------------------------------------------
numpydoc_class_members_toctree = False


def linkcode_resolve(domain, info):
    """Generates the [source] button next to class and function
    definitions in the docs.

    Based on implementation here:
    https://github.com/numpy/numpy/blob/main/doc/source/conf.py

    Parameters
    ----------
    domain : str
        The language domain the object is in. py, c, cpp, or javascript
    info : dict
        Has the following keys guaranteed to be present (dep on domain):
          py: module (name of mod), fullname (name of obj)
          c: names (list of names for the obj)
          cpp: names (list of names for the obj)
          javascript: object (name of obj), fullname (name of the item)

    Returns
    -------
    url : str or None
        The url where the [source] button should link. returns None
        if no link is to be added.
    """

    # LEAVE THE PRINTS! THEY'RE A CRUDE LOGGING/TESTING TOOL!

    # We only want to link python code.
    if domain != 'py':
        print("LJLErr: found a non-python object!")
        print("domain", domain)
        print("info", info)
        return None
    if not info['module']:
        print("LJLErr: found a non-module object!")
        print("info")
        return None

    # Now, finding the correct URL. This depends on 2 separate parts:
    # - finding the correct file (webpage)
    # - finding the relevant line numbers
    # To find these, we first grab the corresponding python object of
    # whatever this function is called for.

    # Import module
    modname = info['module']
    fullname = info['fullname']

    submod = sys.modules.get(modname)
    if submod is None:
        print("LJLErr: submod not found!")
        print("modname", modname)
        # try:
        #     __import__(modname)
        #     submod = sys.modules[modname]
        # except ImportError:
        return None

    # Get the specific function/class object
    obj = submod
    for part in fullname.split('.'):
        try:
            obj = getattr(obj, part)
        except AttributeError:
            return None
        except Exception as ex:
            print("LJLErr: weird exception caught!")
            print(ex)
            return None

    # Unwrap any decorators to find the real source file
    try:
        obj = inspect.unwrap(obj)
    except Exception:
        pass

    # Find file path and line numbers:
    try:
        fn = inspect.getsourcefile(obj)  # get abs path of source file
    except Exception:
        fn = None
    if not fn:
        return None

    module = inspect.getmodule(obj)  # get module
    if module is not None and not module.__name__.startswith("GMAP"):
        print("LJLErr: module from outside GMAP called")
        print("obj", obj)
        return None

    # get module path relative to repo
    fn_path = Path(fn)
    rel_fn = str(fn_path.relative_to(repopath))

    # get line numbers
    try:
        source, lineno = inspect.getsourcelines(obj)
        linespec = f"#L{lineno}-L{lineno + len(source) - 1}"
    except Exception:
        # Fallback for properties or objects where line lookup fails
        linespec = ""
        print("LJLErr: Linespec failed")
        print(fullname)

    # Build GitHub URL
    url = f"https://github.com/kimvana/GMAP_clone/blob/{source_branch}/{rel_fn}{linespec}"
    return url


def release_download_overview():
    subprocess.run(["git", "fetch", "--tags"])
    try:
        tags_raw = subprocess.check_output(
            ["git", "tag", "-l", "--sort=-v:refname"]).decode("utf-8")
        tags = [t.strip() for t in tags_raw.split("\n") if t.strip()]
    except Exception:
        tags = []

    base_url = "https://github.com/Kimvana/GMAP_clone/releases/tag/"
    # return_string = (
    #     "Download link for most recent patch of each version:\n\n"
    # )

    # Generate .rst
    # prev_version = ("x", "x")
    # for tag in tags:
    #     tag_split = tag.split("-")
    #     tag_split = tag_split[0].split(".") + tag_split[1:]
    #     version = (tag_split[0], tag_split[1])
    #     if version != prev_version:
    #         if version[0] != prev_version[0]:
    #             return_string += (
    #                 f"\nGMAP {version[0]}\n{'-' * 30}\n\n"
    #                 ".. hlist::\n    :columns: 4\n\n"
    #             )
    #         return_string += (
    #             f"    - `Version {tag[1:]} <{base_url}{tag}>`__\n"
    #         )
    #     prev_version = version

    # group the tags
    tags_sorted = []
    prev_version = ("x", "x")
    for tag in tags:
        tag_split = tag.split("-")
        tag_split = tag_split[0].split(".") + tag_split[1:]
        version = (tag_split[0], tag_split[1])
        if version != prev_version:
            if version[0] != prev_version[0]:
                tags_sorted.append([version[0], []])
            tags_sorted[-1][1].append(tag)
        prev_version = version

    # build html
    return_string = ""
    for (headertext, subtags) in tags_sorted:
        return_string += (
            f'<section id="gmap-{headertext}">\n'
            f'  <h3>\n'
            f'    "GMAP {headertext}"\n'
            f'    <a class="headerlink" href="#gmap-{headertext}" title="Link to this heading">#</a>\n'
            f'  </h3>\n'
            f'  <table class="hlist">\n'
            f'    <tbody>\n'
            f'      <tr>\n'
        )

        # determine number of items per column
        n_total = len(subtags)
        n_items = n_total // 4
        if n_total % 4 != 0:
            n_items += 1

        for colnum in range(4):
            return_string += (
                f'        <td>\n'
                f'          <ul class="simple">\n'
            )

            for tag in subtags[(n_items * colnum):(n_items * (colnum + 1))]:
                return_string += (
                    f'            <li>\n'
                    f'              ::marker\n'
                    f'              <p>\n'
                    f'                <a class="reference external" href="{base_url}{tag}>Version {tag[1:]}</a>"\n'
                    f'              </p>\n'
                    f'            </li>\n'
                )

            return_string += (
                "          </ul>\n"
                "        </td>\n"
            )

        return_string += (
            "      </tr>\n"
            "    </tbody>\n"
            "  </table>\n"
            "</section>\n"
        )

    return return_string


# linkfile = Path(__file__).parent / "linkfile.rst"
# with open(linkfile, "w", encoding="utf-8") as fhand:
#     fhand.write(release_download_overview())

linkfile = Path(__file__).parent / "_static" / "release_download_overview.html"
with open(linkfile, "w", encoding="utf-8") as fhand:
    fhand.write(release_download_overview())
