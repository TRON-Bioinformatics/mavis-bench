#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# svmerge documentation build configuration file, created by
# sphinx-quickstart on Wed Oct 19 10:47:56 2016.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'
import sys
import os
import sphinx_rtd_theme
import subprocess
import glob
import re

d = os.path.abspath('./../../')
print(d)
sys.path.insert(0, d)
from mavis import __version__
from mavis.config import REFERENCE_DEFAULTS
from mavis.schedule.constants import OPTIONS as SUBMIT_OPTIONS
from mavis.summary.constants import DEFAULTS as SUMMARY_DEFAULTS
from mavis.pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from mavis.validate.constants import DEFAULTS as VALIDATION_DEFAULTS
from mavis.annotate.constants import DEFAULTS as ANNOTATION_DEFAULTS
from mavis.cluster.constants import DEFAULTS as CLUSTER_DEFAULTS
from mavis.illustrate.constants import DEFAULTS as ILLUSTRATION_DEFAULTS
from mavis.util import ENV_VAR_PREFIX

def preprocess():

    dirname = os.path.dirname(os.path.abspath(__file__))

    sys.path.insert(0, os.path.join(dirname, '../..'))

    # auto build the other documentation
    subprocess.check_call('sphinx-apidoc -f -P -o {} {} --separate'.format(os.path.join(dirname, 'auto'), os.path.join(dirname, '../..', 'mavis')), shell=True)

    # now we need to add showing only select special members
    for filename in glob.glob(os.path.join(dirname, 'auto', '*.rst')):
        # now open and read the file
        lines = []
        with open(filename, 'r') as fh:
            lines = fh.readlines()
        with open(filename, 'w') as fh:
            saw_automodule = False
            for line in lines:
                if re.match(r'^\.\.\s+automodule::\s+.*$', line):
                    fh.write(line)
                    fh.write('    :special-members: __and__, __or__, __xor__, __len__, __sub__, __add__\n')
                elif re.match(r'(\S+)\.(\S+)\s+(module|package)', line):
                    m = re.match(r'(\S+)\.(\S+)\s+(module|package)', line)
                    if m.group(1) == 'package':
                        line = re.sub(r'(\S+)\.(\S+)\s+(package)', r'\g<2> package', line)
                    else:
                        line = re.sub(r'(\S+)\.(\S+)\s+(module)', r'\g<2> module', line)
                    fh.write(line)
                else:
                    fh.write(line)

    fname = os.path.join(dirname, 'config_settings_glossary.rst')
    print('writing:', fname)
    with open(fname, 'w') as fh:
        fh.write('Configurable Settings\n')
        fh.write('+' * 50)
        tab = ' ' * 4
        fh.write('\n\n.. glossary::\n{}:sorted:\n\n'.format(tab))
        glossary = {}
        CUSTOM_TYPES = {
            'cast_boolean': 'bool',
            'float_fraction': '~mavis.constants.float_fraction',
            'ChrListString': '~mavis.util.ChrListString'
        }
        for namespace in [
            SUBMIT_OPTIONS,
            REFERENCE_DEFAULTS,
            SUMMARY_DEFAULTS,
            PAIRING_DEFAULTS,
            ANNOTATION_DEFAULTS,
            VALIDATION_DEFAULTS,
            CLUSTER_DEFAULTS,
            ILLUSTRATION_DEFAULTS
        ]:
            for term, value in namespace.items():
                typ = namespace.type(term).__name__
                typ = CUSTOM_TYPES.get(typ, typ)
                defn = ':class:`{}` - {}. The corresponding environment variable is ``{}{}`` and the default value is ``{}``'.format(
                    typ, re.sub(r'\.?$', '', namespace.define(term, '')).capitalize(), ENV_VAR_PREFIX, term.upper(), repr(value))
                try:
                    defn += '. Accepted values include: {}'.format(', '.join(['``{}``'.format(repr(v)) for v in namespace.type(term).values()]))
                except AttributeError:
                    pass
                glossary[term] = defn
        for term, defn in sorted(glossary.items()):
            fh.write('\n{}{}\n'.format(tab, term))
            fh.write('{}{}{}\n'.format(tab, tab, defn))

preprocess()

d = os.path.abspath('../../bin')
sys.path.insert(0, d)

print('added to path:', d)

autoclass_content = 'both'
autodoc_default_flags = ['show-inheritance', 'members']
todo_include_todos = True
# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx.ext.intersphinx',
    'm2r'
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3.6', None),
    'pysam': ('http://pysam.readthedocs.io/en/latest/', None),
    'networkx': ('https://networkx.readthedocs.io/en/stable/', None),
    'svgwrite': ('https://pythonhosted.org/svgwrite/', None)
}

mathjax_path="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"

html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
max_depth = -1

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
numfig = True
numfig_format = {'figure': 'Figure %s.', 'table': 'Table: %s', 'code-block': ''}
source_suffix = ['.rst', '.md']

# The encoding of source files.
#
# source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = 'MAVIS'
copyright = '2017, creisle'
author = 'creisle'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = __version__
# The full version, including alpha/beta/rc tags.
#release = '1.0.0'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = []


# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []

# If true, keep warnings as "system message" paragraphs in the built documents.
# keep_warnings = False

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

# The name of an image file (relative to this directory) to use as a favicon of
# the docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#
# html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Output file base name for HTML help builder.
htmlhelp_basename = 'MAVISdoc'


def setup(app):
    app.add_css_file('custom.css')
