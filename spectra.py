import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
import numpy as np
from astroquery.simbad import Simbad
import time

QUERY_DELAY = 0.6
last_query_time = 0

spectrum_mode = "absorption"
autocomplete_cache = {}

current_lines = []
current_text = ""
# -------------------
# Data
# -------------------
spectral_data = {
    "H": [656, 486, 434, 410],
    "He": [501, 587, 706, 728],
    "Li": [670, 610, 497],
    "Be": [313, 234],
    "B": [249, 208],
    "C": [156, 165, 247],
    "N": [149, 174, 178],
    "O": [557, 630, 636],
    "F": [640, 650],
    "Fe": [438, 527, 533],
    "Ne": [540, 585, 703, 743],
    "Na": [589, 590],
    "Mg": [285, 383, 518],
    "Al": [396, 394, 309],
    "Si": [288, 250],
    "P": [177, 179],
    "S": [182, 186],
    "Cl": [519, 550],
    "Ar": [750, 842],
    "K": [766, 770],
    "Ca": [422, 393],
}

element_names = {
    "H": "Hydrogen", "He": "Helium", "Li": "Lithium", "Be": "Beryllium",
    "B": "Boron", "C": "Carbon", "N": "Nitrogen", "O": "Oxygen",
    "F": "Fluorine", "Ne": "Neon", "Na": "Sodium", "Mg": "Magnesium",
    "Al": "Aluminum", "Si": "Silicon", "P": "Phosphorus", "S": "Sulfur",
    "Cl": "Chlorine", "Ar": "Argon", "K": "Potassium", "Ca": "Calcium",
    "Fe": "Iron"
}

astro_objects = {
    "Sun": ["H", "He", "O", "C", "Fe", "Ca"],
    "Sirius":["H", "He"],
    "Betelgeuse":["H", "He", "C", "O", "Si"],
    "Earth": ["O", "N", "Ar", "C"],
    "Mars": ["O", "C", "Fe"],
    "Jupiter":["H", "He"],
    "Saturn":["H", "He"],
}

stellar_temperatures = {
    "O": "30,000+ K",
    "B": "10,000–30,000 K",
    "A": "7,500–10,000 K",
    "F": "6,000–7,500 K",
    "G": "5,200–6,000 K",
    "K": "3,700–5,200 K",
    "M": "2,400–3,700 K"
}

stellar_classes = {
    "O": ["He"], "B": ["He", "H"], "A": ["H"],
    "F": ["H", "Ca", "Fe"], "G": ["H", "Fe", "Ca"],
    "K": ["Fe", "Ca"], "M": ["O", "C"]
}

# -------------------
# SIMBAD API
# -------------------
def search_object(name):

    simbad = Simbad()
    simbad.add_votable_fields("otype","sp_type","plx_value","ra","dec")

    try:
        result = simbad.query_object(name)

        if result is None:
            return None

        return {
            "type": str(result["otype"][0]),
            "spectral_type": str(result["sp_type"][0]),
            "parallax": str(result["plx_value"][0]),
            "ra": str(result["ra"][0]),
            "dec": str(result["dec"][0])
        }

    except Exception as e:
        print("SIMBAD error:", e)
        return None
# -------------------
# Spectrum functions
# -------------------
def wavelength_to_rgb(wl, gamma=0.8):

    if wl < 380 or wl > 750:
        return (0,0,0)

    if wl < 440:
        attenuation=0.3 + 0.7*(wl-380)/(440-380)
        R = (-(wl-440)/(440-380)*attenuation)**gamma
        G = 0
        B = attenuation**gamma

    elif wl < 490:
        R = 0
        G = ((wl-440)/(490-440))**gamma
        B = 1

    elif wl < 510:
        R = 0
        G = 1
        B = (-(wl-510)/(510-490))**gamma

    elif wl < 580:
        R = ((wl-510)/(580-510))**gamma
        G = 1
        B = 0

    elif wl < 645:
        R = 1
        G = (-(wl-645)/(645-580))**gamma
        B = 0

    else:
        attenuation = 0.3 + 0.7*(750-wl)/(750-645)
        R = attenuation**gamma
        G = 0
        B = 0

    return (R,G,B)

def draw_spectrum_bg():

    wl = np.linspace(380,750,1000)
    gradient = np.zeros((1,len(wl),3))

    for i,w in enumerate(wl):
        gradient[0,i] = wavelength_to_rgb(w)

    ax.imshow(gradient,aspect="auto",extent=[380,750,0,1])


def parallax_to_distance(plx):

    try:
        plx = float(plx)
        if plx > 0:
            return round(1000/plx,2)
    except:
        pass

    return "Unknown"


def plot_spectrum(elements, title):

    ax.clear()
    draw_spectrum_bg()

    lines = []

    spectral_lines = sorted(
        (wl, el)
        for el in elements
        if el in spectral_data
        for wl in spectral_data[el]
        if 380 <= wl <= 750
    )

    for wl, el in spectral_lines:

        color = "black"
        if spectrum_mode == "emission":
            color = wavelength_to_rgb(wl)

        line = ax.axvline(wl, color=color, linewidth=2, picker=5)

        line.element = el
        line.wavelength = wl

        lines.append(line)

    ax.set_xlim(380,750)
    ax.set_ylim(0,1)
    ax.set_title(title)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_yticks([])

    fig.canvas.draw_idle()

    return lines

# -------------------
# GUI
# -------------------
fig, ax = plt.subplots(figsize=(10,4))
plt.subplots_adjust(top=0.8,bottom=0.2)

info_text = fig.text(
    0.5,0.05,"",
    ha="center",
    fontsize=10,
    bbox=dict(boxstyle="round",fc="white")
)

# Hover annotation
hover_annotation = ax.annotate(
    "",
    xy=(0,0),
    xytext=(15,15),
    textcoords="offset points",
    bbox=dict(boxstyle="round,pad=0.3", fc="black", alpha=0.8),
    color="white",
    fontsize=11
)

hover_annotation.set_visible(False)

# -------------------
# Events
# -------------------
def autocomplete(prefix):

    if prefix in autocomplete_cache:
        return autocomplete_cache[prefix]

    try:
        simbad = Simbad()
        result = simbad.query_objectids(prefix)

        if result is None:
            return []

        names = [
            str(r["ID"]).replace("NAME ","")
            for r in result[:5]
        ]

        autocomplete_cache[prefix] = names

        return names

    except:
        return []

def on_click(event):

    if event.inaxes != ax or event.xdata is None:
        return

    threshold = 2

    for line in current_lines:

        if abs(event.xdata - line.wavelength) < threshold:

            name = element_names.get(line.element, line.element)

            print(f"Selected {name} at {line.wavelength} nm")

def on_hover(event):

    if event.inaxes != ax or event.xdata is None:
        hover_annotation.set_visible(False)
        fig.canvas.draw_idle()
        return

    closest_line = None
    closest_distance = float("inf")

    for line in current_lines:
        distance = abs(event.xdata - line.wavelength)

        if distance < closest_distance:
            closest_distance = distance
            closest_line = line

    HOVER_THRESHOLD = 2

    # reset all line widths first
    for l in current_lines:
        if l is not None:
            l.set_linewidth(2)

    if closest_line is not None and closest_distance < HOVER_THRESHOLD:
        name = element_names.get(closest_line.element, closest_line.element)

        hover_annotation.xy = (closest_line.wavelength, 0.9)
        hover_annotation.set_position((15, 15))
        hover_annotation.set_text(f"{name}\n{closest_line.wavelength} nm")

        hover_annotation.set_visible(True)

        closest_line.set_linewidth(4)

    else:
        hover_annotation.set_visible(False)

    fig.canvas.draw_idle()


'''# TextBox
text_box_ax = plt.axes([0.35,0.9,0.4,0.05])
text_box = TextBox(text_box_ax,'Enter Planet/Star or Class: ')
suggestion_text = fig.text(
    0.5,
    0.82,
    "",
    fontsize=9,
    color="blue"
)
text_box.on_text_change(on_text_change)'''

def submit(text):

    global current_lines

    text = text.strip()

    if text.upper() in stellar_classes:
        elements = stellar_classes[text.upper()]
        current_lines[:] = plot_spectrum(elements,f"{text.upper()}-Type Star")
        info_text.set_text(f"Local stellar class: {text.upper()}")
        return

    if text in astro_objects:
        elements = astro_objects[text]
        current_lines[:] = plot_spectrum(elements,f"{text} Spectrum")
        info_text.set_text(f"Local object: {text}")
        return

    result = search_object(text)

    if not result:
        info_text.set_text("Unknown object")
        return

    spectral_type = result.get("spectral_type","Unknown")

    star_class = spectral_type[0] if spectral_type != "Unknown" else "G"

    elements = stellar_classes.get(star_class, stellar_classes["G"])

    current_lines[:] = plot_spectrum(elements, f"{text} ({star_class}-Type)")

    distance = parallax_to_distance(result.get("parallax"))

    info_text.set_text(
        f"Type: {result.get('type')} | "
        f"Spectral: {spectral_type} | "
        f"Distance: {distance} pc"
    )

fig.canvas.mpl_connect("motion_notify_event", on_hover)
fig.canvas.mpl_connect("button_press_event", on_click)

text_box_ax = plt.axes([0.35,0.9,0.4,0.05])
text_box = TextBox(text_box_ax,"Enter Planet/Star or Class: ")
text_box.on_submit(submit)

plt.show()