multishift stack -i al.vasp al.vasp al.vasp al.vasp al.vasp al.vasp al.vasp al.vasp -o slab.vasp
multishift chain --input slab.vasp --output chain --cleave -0.15 -0.1 -0.05 0.0 0.05 0.1 0.15 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.25 1.5 1.75 2.0 2.5 3.0 3.5 4.0 6.0 8.0 --shift 12 12
