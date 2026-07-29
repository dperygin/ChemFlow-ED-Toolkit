from pymol import cmd
import math
import csv

def distance(a, b):
    return math.sqrt(sum([(a[i]-b[i])**2 for i in range(3)]))

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

def is_hbond_donor(atom):
    return atom["elem"] in ["N","O"]

def is_hbond_acceptor(atom):
    return atom["elem"] in ["N","O"]

def is_hydrophobic(atom):
    return atom["elem"] == "C"

def is_aromatic(atom):
    # crude approximation (PyMOL lacks direct aromatic flag here)
    return atom["elem"] == "C"

def find_interactions(lig_sel, rec_sel, cutoff=5.0, outfile="interactions.csv"):

    lig_atoms = get_atoms(lig_sel)
    rec_atoms = get_atoms(rec_sel)

    interactions = []

    for la in lig_atoms:
        for ra in rec_atoms:
            d = distance(la["coord"], ra["coord"])
            if d > cutoff:
                continue

            interaction_type = None

            # H-bond donor (ligand donor)
            if is_hbond_donor(la) and is_hbond_acceptor(ra) and d < 3.5:
                interaction_type = "Hbond_donor"

            # H-bond acceptor (ligand acceptor)
            elif is_hbond_acceptor(la) and is_hbond_donor(ra) and d < 3.5:
                interaction_type = "Hbond_acceptor"

            # Ionic
            elif la["charge"] != 0 and ra["charge"] != 0 and (la["charge"] * ra["charge"] < 0) and d < 5.0:
                interaction_type = "ionic"

            # Hydrophobic
            elif is_hydrophobic(la) and is_hydrophobic(ra) and d < 4.5:
                interaction_type = "hydrophobic"

            # Pi-pi (very simplified)
            elif is_aromatic(la) and is_aromatic(ra) and d < 5.0:
                interaction_type = "pi-pi"

            if interaction_type:
                interactions.append({
                    "lig_atom": f"{la['resn']}{la['resi']}/{la['name']}/{la['index']}",
                    "rec_atom": f"{ra['resn']}{ra['resi']}/{ra['name']}/{ra['index']}",
                    "type": interaction_type,
                    "dist": d
                })

    # keep shortest distance per ligand atom + interaction type
    best = {}
    for i in interactions:
        key = (i["lig_atom"], i["type"])
        if key not in best or i["dist"] < best[key]["dist"]:
            best[key] = i

    final = list(best.values())

    # print to console
    print("LigandAtom,ReceptorAtom,InteractionType,Distance")
    for i in final:
        print(f"{i['lig_atom']},{i['rec_atom']},{i['type']},{i['dist']:.2f}")

    # write CSV
    with open(outfile, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["LigandAtom","ReceptorAtom","InteractionType","Distance"])
        for i in final:
            writer.writerow([i["lig_atom"], i["rec_atom"], i["type"], f"{i['dist']:.3f}"])

    print(f"\nSaved to {outfile}")

# run it
find_interactions("docked_ligand_0001 and sele", "receptor and sele")
