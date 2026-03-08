import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
import numpy as np
from astroquery.simbad import Simbad
import time

last_query_time = 0
QUERY_DELAY = 0.6

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
spectrum_mode = "absorption"
autocomplete_cache = {}
# -------------------
# SIMBAD API
# -------------------
def search_object(name):
    try:
        simbad = Simbad()
        simbad.add_votable_fields('otype','sp_type','plx_value', 'ra', 'dec')

        result = simbad.query_object(name)

        if result is None:
            print("Object not found.")
            return None

        obj_type = result['otype'][0] if 'otype' in result.colnames else "Unknown"
        spectral_type = result['sp_type'][0] if 'sp_type' in result.colnames else "Unknown"
        parallax = result['plx_value'][0] if 'plx_value' in result.colnames else "Unknown"
        ra = result['ra'][0] if 'ra' in result.colnames else "Unknown"
        dec = result['dec'][0] if 'dec' in result.colnames else "Unknown"

        print(result)
        print(result.colnames)

        return {
            "type": str(obj_type),
            "spectral_type": str(spectral_type),
            "parallax": str(parallax),
            "ra": str(ra),
            "dec": str(dec)
        }

    except Exception as e:
        print("Search error:", e)
        return None

# -------------------
# Spectrum functions
# -------------------
def wavelength_to_rgb(wl, gamma=0.8):
    if wl < 380 or wl > 750:
        return 0,0,0
    if wl<440:
        att=0.3+0.7*(wl-380)/(440-380)
        R=(-(wl-440)/(440-380)*att)**gamma; G=0; B=att**gamma
    elif wl<490: R=0; G=((wl-440)/(490-440))**gamma; B=1
    elif wl<510: R=0; G=1; B=(-(wl-510)/(510-490))**gamma
    elif wl<580: R=((wl-510)/(580-510))**gamma; G=1; B=0
    elif wl<645: R=1; G=(-(wl-645)/(645-580))**gamma; B=0
    else: att=0.3+0.7*(750-wl)/(750-645); R=att**gamma; G=0; B=0
    return R,G,B

def draw_spectrum_bg():
    wl = np.linspace(380,750,1000)
    grad = np.zeros((1,len(wl),3))
    for i,w in enumerate(wl): grad[0,i]=wavelength_to_rgb(w)
    ax.imshow(grad,aspect='auto',extent=[380,750,0,1])

def parallax_to_distance(plx):
    try:
        plx = float(plx)
        if plx > 0:
            return round(1000/plx,2)
    except:
        return "Unknown"

def plot_absorption(elements, title):
    ax.clear()
    draw_spectrum_bg()

    lines = []

    for el in elements:
        if el in spectral_data:
            for wl in spectral_data[el]:

                swl = wl

                if 380 <= swl <= 750:

                    color = 'black'
                    if spectrum_mode == "emission":
                        color = wavelength_to_rgb(swl)

                    line = ax.axvline(swl, color=color, linewidth=2, picker=5)

                    line.element = el
                    line.wavelength = swl
                    lines.append(line)

                    ax.text(
                        2 + swl,
                        0.5,
                        el,
                        rotation=45,
                        ha='center',
                        va='bottom',
                        fontsize=15
                    )

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
info_text = fig.text(
    .5, 0.05,
    "",
    ha="center",
    fontsize=10,
    bbox=dict(boxstyle="round", fc="white")
)
plt.subplots_adjust(top=0.8,bottom=0.2)

# Hover annotation
hover_annotation = ax.annotate("", xy=(0,0), xytext=(10,10), textcoords='offset points',
                               bbox=dict(boxstyle='round',fc='white'))
hover_annotation.set_visible(False)

# -------------------
# Events
# -------------------

def on_text_change(text):
    global last_query_time

    text = text.strip()

    # Don't search until at least 3 characters
    if len(text) < 3:
        suggestion_text.set_text("")
        fig.canvas.draw_idle()
        return

    now = time.time()

    if now - last_query_time < QUERY_DELAY:
        return

    last_query_time = now

    suggestions = autocomplete(text)

    if suggestions:
        suggestion_text.set_text("Suggestions: " + ", ".join(suggestions))
    else:
        suggestion_text.set_text("")

    fig.canvas.draw_idle()

def autocomplete(prefix):

    if prefix in autocomplete_cache:
        return autocomplete_cache[prefix]

    try:
        simbad = Simbad()
        result = simbad.query_objectids(prefix)

        if result is None:
            return []

        names = []

        for r in result[:5]:
            name = str(r["ID"]).replace("NAME ", "")
            names.append(name)

        autocomplete_cache[prefix] = names

        return names

    except Exception:
        return []


current_lines = []
current_text = ""

def on_hover(event):
    if event.inaxes != ax or event.xdata is None:
        hover_annotation.set_visible(False)
        fig.canvas.draw_idle()
        return
    threshold = 2
    for line in current_lines:
        if abs(event.xdata - line.wavelength) < threshold:
            full_name = element_names.get(line.element,line.element)
            hover_annotation.xy = (line.wavelength,1)
            hover_annotation.set_text(full_name)
            hover_annotation.set_visible(True)
            fig.canvas.draw_idle()
            return
    hover_annotation.set_visible(False)
    fig.canvas.draw_idle()

fig.canvas.mpl_connect("motion_notify_event", on_hover)

# TextBox
text_box_ax = plt.axes([0.35,0.9,0.4,0.05])
text_box = TextBox(text_box_ax,'Enter Planet/Star or Class: ')
suggestion_text = fig.text(
    0.5,
    0.82,
    "",
    fontsize=9,
    color="blue"
)
text_box.on_text_change(on_text_change)

def submit(text):
    global current_lines,current_text
    current_text = text.strip()

    # Local dictionary
    if current_text.upper() in stellar_classes:
        elements = stellar_classes[current_text.upper()]
        current_lines[:] = plot_absorption(elements,f"{current_text.upper()}-Type Star Spectrum")
        info_text.set_text(f"Local stellar class: {current_text.upper()}")
        fig.canvas.draw_idle()
        return
    if current_text in astro_objects:
        elements = astro_objects[current_text]
        current_lines[:] = plot_absorption(elements,f"{current_text} Spectrum")
        info_text.set_text(f"Local object: {current_text}")
        fig.canvas.draw_idle()
        return
    # SIMBAD search
    result = search_object(current_text)

    if result:
        print("Found:", result)

        spectral_type = result.get("spectral_type", "Unknown")

        temp = "Unknown"
        distance = "Unknown"

        if result.get("parallax") != "Unknown":
            distance = parallax_to_distance(result["parallax"])

        if spectral_type != "Unknown" and len(spectral_type) > 0:
            star_class = spectral_type[0].upper()

            temp = stellar_temperatures.get(star_class, "Unknown")

            if star_class in stellar_classes:
                elements = stellar_classes[star_class]
                current_lines[:] = plot_absorption(elements, f"{current_text} ({star_class}-Type)")
        else:
            current_lines[:] = plot_absorption(stellar_classes["G"], f"{current_text} (Unknown Star)")

        info_text.set_text(
            f"Type: {result.get('type', 'Unknown')} | "
            f"Spectral: {result.get('spectral_type', 'Unknown')} | "
            f"Temp: {temp} | "
            f"Parallax: {result.get('parallax', 'Unknown')} | "
            f"Distance: {distance} pc"
        )
        fig.canvas.draw_idle()
    else:
        print("Unknown object.")
        info_text.set_text("Unknown object")
        fig.canvas.draw_idle()

text_box.on_submit(submit)

plt.show()