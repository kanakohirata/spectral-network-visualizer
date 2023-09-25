import base64
from io import BytesIO
from logging import getLogger
import re
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import GetPeriodicTable as pt
from rdkit.Chem.Draw import rdMolDraw2D

logger = getLogger(__name__)


def generate_inchikey(structure):
    """Generate InChIKey using RDKit

    Parameters
    ----------
    structure : str
        InChI or SMILES.

    Returns
    -------
    str
        InChIKey or empty string.

    """

    if structure.startswith('InChI'):
        mol = Chem.MolFromInchi(structure)
    else:
        mol = Chem.MolFromSmiles(structure)

    if mol:
        inchikey = Chem.MolToInchiKey(mol)
        return inchikey

    return ''


def structure_matching(structure1, structure2, image_width, image_height):
    """

    Parameters
    ----------
    structure1 : str
        InChI or SMILES.

    structure2 : str
        InChI or SMILES.

    image_width : int
        Width (px) of SVG image.

    image_height : int
        Height (px) of SVG image.

    Returns
    -------
    list
        List of SVG(str) for structure1 with highlighted substructure that match structure2.

    """

    if structure1.startswith('InChI'):
        mol1 = Chem.MolFromInchi(structure1)
    else:
        mol1 = Chem.MolFromSmiles(structure1)

    if structure2.startswith('InChI'):
        mol2 = Chem.MolFromInchi(structure2)
    else:
        mol2 = Chem.MolFromSmiles(structure2)

    if mol1 and mol2:
        matches = mol1.GetSubstructMatches(mol2)

        l_svg = []
        for match in matches:
            # Generate container
            drawer = rdMolDraw2D.MolDraw2DSVG(image_width, image_height)

            # Set options
            options = drawer.drawOptions()
            options.clearBackground = False

            # Pass mol to container and get svg text
            drawer.DrawMolecule(mol1, highlightAtoms=match)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()

            # Edit header in svg
            size_inline = f"width='{image_width}px' height='{image_height}px'"
            svg = svg.replace(size_inline, "")
            svg_header, svg_body = svg.split("<!-- END OF HEADER -->\n")
            l_svg_header_line = [line.strip() for line in svg_header.split("\n")[1:]]
            new_svg_header = "\n".join(l_svg_header_line)
            svg = new_svg_header + svg_body

            l_svg.append(svg)

        return l_svg


def get_mol_structure_2dsvg(mol, width, height) -> str:
    """SVG string of 2D structure.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol object.

    width : int
        Width (px) of SVG image.

    height : int
        Height (px) of SVG image.

    Returns
    -------
    str
        SVG of 2D structure for the mol object.

    """
    # Generate container
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)

    # Set options
    options = drawer.drawOptions()
    options.clearBackground = False

    # Pass mol to container and get svg text
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # Edit header in svg
    size_inline = f"width='{width}px' height='{height}px'"
    svg = svg.replace(size_inline, "")
    svg_header, svg_body = svg.split("<!-- END OF HEADER -->\n")
    l_svg_header_line = [line.strip() for line in svg_header.split("\n")[1:]]
    new_svg_header = "\n".join(l_svg_header_line)
    svg = new_svg_header + svg_body

    return svg


def get_mol_structure_2dsvg_base64(mol, width, height) -> str:
    """SVG string of 2D structure.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol object.

    width : int
        Width (px) of SVG image.

    height : int
        Height (px) of SVG image.

    Returns
    -------
    str
        Base64 encoded SVG of 2D structure for the mol object.

    """
    # Generate container
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)

    # Set options
    options = drawer.drawOptions()
    options.clearBackground = False

    # Pass mol to container and get svg text
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # Edit header in svg
    size_inline = f"width='{width}px' height='{height}px'"
    svg = svg.replace(size_inline, "")
    svg_header, svg_body = svg.split("<!-- END OF HEADER -->\n")
    l_svg_header_line = [line.strip() for line in svg_header.split("\n")[1:]]
    new_svg_header = "\n".join(l_svg_header_line)
    svg = new_svg_header + svg_body

    buffered = BytesIO()
    buffered.write(str.encode(svg))
    _svg_base64 = base64.b64encode(buffered.getvalue())
    svg_base64 = f"data:image/svg+xml;base64,{repr(_svg_base64)[2:-1]}"

    return svg_base64


def parse_formula(formula):
    """Separate formula into components

    Parameters
    ----------
    formula : str
        Molecular formula e.g. C9H17NO4

    Returns
    -------
    set
        Set of elements e.g. {C9, H17, N, O4}
    """

    count = 0
    s_element = set()

    while formula:
        count += 1
        m = re.match(r'[A-Z][a-z]?\d*|\((?:[^()]*(?:\(.*\))?[^()]*)+\)\d+', formula)
        if not m:
            break
        s_element.add(m.group())
        formula = formula.lstrip(m.group())

        if count == 120:
            break

    return s_element


def calc_ion_mass(ion_form):
    """Calculate exact mass form ion formula

    Parameters
    ----------
    ion_form : str
        Charged form of ion e.g. [C19H28O-H]+

    Returns
    -------
    float
        Theoretical ion mass
    """

    emass = 0.00054857990924
    charge = ion_form.split(']')[1]
    ion_form = ion_form.split('[', 1)[1].split(']')[0]

    formula = re.split(r'[+-]', ion_form, 1)[0]
    formula_diff = ion_form.replace(formula, '', 1)

    logger.debug(f'formula: {formula}, formula_diff: {formula_diff}')

    d_elem_vs_num = {}

    elems = re.findall(r'[A-Z][a-z]?\d*', formula)

    if not elems:
        raise ValueError(f'{formula} is incorrect formula.')

    for elem in elems:
        symb = re.search(r'[A-Z][a-z]?', elem).group()
        num = 1
        if re.search(r'\d+', elem):
            num = int(re.search(r'\d+', elem).group())

        d_elem_vs_num[symb] = num

    elems_diff = re.findall(r'[+-]\d*[A-Z][a-z]?', formula_diff)
    for elem_diff in elems_diff:
        symb = re.search(r'[A-Z][a-z]?', elem_diff).group()
        num = 1
        if re.search(r'\d+', elem_diff):
            num = int(re.search(r'\d+', elem_diff).group())
        if elem_diff.startswith('-'):
            num *= -1

        if symb in d_elem_vs_num.keys():
            d_elem_vs_num[symb] += num
        else:
            d_elem_vs_num[symb] = num

    # logger.debug(f'd_elem_vs_num: {d_elem_vs_num}')

    ion_mass = 0
    for symb, num in d_elem_vs_num.items():
        elem_mass = pt().GetMostCommonIsotopeMass(symb) * num
        # logger.debug(f'{symb} * {num} = {elem_mass}')
        ion_mass += elem_mass

    charge_num = 1
    if re.search(r'\d+', charge):
        charge_num = int(re.search(r'\d+', elem).group())

    if charge.endswith('+'):
        ion_mass -= emass * charge_num
    elif charge.endswith('-'):
        ion_mass += emass * charge_num
    else:
        logger.warning(f'{ion_form} is not charged. Calculation of electrons mass is not performed.')

    # logger.debug(ion_mass)

    return ion_mass
