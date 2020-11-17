SHAPE=100
RADIUS=10
VALUE=1

# BEGIN: NO CHANGE

VOLDIR=${ROOTDIR}/vol
DATADIR=${ROOTDIR}/data
IMGDIR=${ROOTDIR}/images
DOCDIR=${ROOTDIR}/doc

FORCEFIELD=amber

PDB2PQR=~/Software/pdb2pqr-1.7/pdb2pqr.py --ff ${FORCEFIELD}

vpath %.ply ${DATADIR}
vpath %.pdb ${DATADIR}
vpath %.pqr ${DATADIR}

# END: NO CHANGE
