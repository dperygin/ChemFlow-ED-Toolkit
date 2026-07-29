from pymol import cmd
import math
import csv

# ---------- BASIC UTILITIES ----------

def distance(a, b):
    return math.sqrt(sum([(a[i]-b[i])**2 for i in range(3)]))

def centroid(coords):
    return [sum(x)/len(coords) for x in zip(*coords)]

# ---------- ATOM EXTRACTION ----------

def get_atoms(selection):
    model = cmd.get_model(selection)
    atoms = []
    for a in model.atom:
        atoms.append({
            "name": a.name,
            "resn": a.resn,
            "resi": a.resi,
            "coord": a.coord,
            "elem": a.symbol,
            "charge": a.formal_charge,
            "index": a.index
        })
    return atoms

# ---------- INTERACTION FLAGS ----------

def is_hbond_donor(atom):
    return atom["elem"] in ["N", "O"]

def is_hbond_acceptor(atom):
    return atom["elem"] in ["N", "O"]

def is_hydrophobic(atom):
    return atom["elem"] == "C"

# ---------- PI-PI DETECTION (ROBUST) ----------

def get_protein_rings(selection):
    aromatic_residues = ["PHE", "TYR", "TRP", "HIS"]
    model = cmd.get_model(selection)
    rings = {}

    for a in model.atom:
        if a.resn in aromatic_residues:
            key = (a.resi, a.resn)
            rings.setdefault(key, []).append(a.coord)

    return rings


def get_ligand_rings(selection):
    model = cmd.get_model(selection)
    coords = [a.coord for a in model.atom if a.symbol == "C"]

    if len(coords) >= 5:
        return {("LIG", "UNK"): coords}
    return {}


def detect_pi_pi(lig_sel, rec_sel):
    interactions = []

    lig_rings = get_ligand_rings(lig_sel)
    rec_rings = get_protein_rings(rec_sel)

    for lkey, lcoords in lig_rings.items():
        lc = centroid(lcoords)

        for rkey, rcoords in rec_rings.items():
            rc = centroid(rcoords)

            d = distance(lc, rc)

            if d < 5.5:
                interactions.append({
                    "lig_atom": f"{lkey[1]}{lkey[0]}/ring",
                    "rec_atom": f"{rkey[1]}{rkey[0]}/ring",
                    "type": "pi-pi",
                    "dist": d,
                    "centroid_pair": (lc, rc)
                })

    return interactions

# ---------- MAIN FUNCTION ----------

def find_interactions(lig_sel, rec_sel, outfile="interactions.csv"):

    lig_atoms = get_atoms(lig_sel)
    rec_atoms = get_atoms(rec_sel)

    interactions = []

    for la in lig_atoms:
        for ra in rec_atoms:

            d = distance(la["coord"], ra["coord"])
            if d > 5.0:
                continue

            interaction_type = None

            if is_hbond_donor(la) and is_hbond_acceptor(ra) and d < 3.5:
                interaction_type = "Hbond_donor"

            elif is_hbond_acceptor(la) and is_hbond_donor(ra) and d < 3.5:
                interaction_type = "Hbond_acceptor"

            elif la["charge"] != 0 and ra["charge"] != 0 and (la["charge"] * ra["charge"] < 0):
                interaction_type = "ionic"

            elif is_hydrophobic(la) and is_hydrophobic(ra) and d < 4.5:
                interaction_type = "hydrophobic"

            if interaction_type:
                interactions.append({
                    "lig_atom": f"{la['resn']}{la['resi']}/{la['name']}/{la['index']}",
                    "rec_atom": f"{ra['resn']}{ra['resi']}/{ra['name']}/{ra['index']}",
                    "type": interaction_type,
                    "dist": d
                })

    interactions.extend(detect_pi_pi(lig_sel, rec_sel))

    best = {}
    for i in interactions:
        key = (i["lig_atom"], i["type"])
        if key not in best or i["dist"] < best[key]["dist"]:
            best[key] = i

    final = list(best.values())

    with open(outfile, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["LigandAtom", "ReceptorAtom", "InteractionType", "Distance"])
        for i in final:
            writer.writerow([i["lig_atom"], i["rec_atom"], i["type"], f"{i['dist']:.3f}"])

    print(f"Saved CSV → {outfile}")

    cmd.delete("int_*")

    for i, interaction in enumerate(final):

        itype = interaction["type"]
        dist = interaction["dist"]

        obj_name = f"int_{i}_{itype}"

        if itype == "pi-pi" and "centroid_pair" in interaction:
            lc, rc = interaction["centroid_pair"]

            pseudo1 = f"pseudo_lig_{i}"
            pseudo2 = f"pseudo_rec_{i}"

            cmd.pseudoatom(pseudo1, pos=lc)
            cmd.pseudoatom(pseudo2, pos=rc)

            cmd.distance(obj_name, pseudo1, pseudo2)

        else:
            lig_index = interaction["lig_atom"].split("/")[-1]
            rec_index = interaction["rec_atom"].split("/")[-1]

            cmd.distance(
                obj_name,
                f"id {lig_index} and docked_ligand_0001",
                f"id {rec_index} and receptor"
            )

        cmd.label(obj_name, f'"{itype} ({dist:.2f} Å)"')

        if "Hbond" in itype:
            cmd.set("dash_color", "cyan", obj_name)
        elif itype == "hydrophobic":
            cmd.set("dash_color", "yellow", obj_name)
        elif itype == "pi-pi":
            cmd.set("dash_color", "magenta", obj_name)
        elif itype == "ionic":
            cmd.set("dash_color", "red", obj_name)

    cmd.group("interactions", "int_*")

    cmd.set("dash_width", 2.5)
    cmd.set("label_size", 16)
    cmd.bg_color("white")

    cmd.show("sticks", lig_sel)
    cmd.show("cartoon", rec_sel)
    cmd.color("cyan", lig_sel)

    cmd.save("interactions_session.pse")
    print("Saved PyMOL session → interactions_session.pse")


find_interactions("docked_ligand_0001", "receptor")
