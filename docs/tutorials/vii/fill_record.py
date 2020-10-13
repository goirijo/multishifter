import json
import os

def load_record(record_file):
    with open(record_file) as f:
        records = json.load(f)

    return records

def write_json(data,target_file):
    with open(target_file,'w') as f:
        json.dump(data,f,indent=4)

    return

def read_oszicar_energy(oszicar_path):
    with open(oszicar_path) as f:
        content = f.readlines()

    return float(content[-1].split()[2])

def fill_record_with_dft(record,use_symmetry=True):
    equivalents = record["equivalents"]

    if use_symmetry:
        for orbit in equivalents:
            prototype = orbit[0]
            protoszi=os.path.join(record["ids"][prototype]["directory"],"OSZICAR")
            energy = read_oszicar_energy(protoszi)

            for i in orbit:
                record["ids"][i]["dft_energy"]=energy

    else:
        for i in record["ids"]:
            oszipath=os.path.join(record["ids"][i]["directory"],"OSZICAR")
            record["ids"][i]["dft_energy"]=read_oszicar_energy(oszipath)

    return record

def main():
    record=load_record("./record.json")
    record=fill_record_with_dft(record)
    write_json(record,"./modified_record.json")

    return

if __name__ == "__main__":
    main()
