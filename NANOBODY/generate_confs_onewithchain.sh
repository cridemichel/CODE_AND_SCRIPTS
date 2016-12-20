SPHERE_DIAM="1.0"
PATCHY_SPHERE_DIAM="0.119"
ANTIGENE_DIAM="0.790"
NANOLEN="2.0"
NANODIAM="1.0"
NUMSPH="10"
NANOPERMPATCH="0.119"
#./create_one_NANOBODY_attached startOneAttch.cnf 1 1 1 0 0 0 -1 10 2.0 2.0 0.8 0.45 1.0 $SPHERE_DIAM $PATCHY_SPHERE_DIAM $ANTIGENE_DIAM 0
./create_one_NANOBODY_attached startOneAttch.cnf 1 1 1 0 0 0 -1 $NUMSPH $NANODIAM $NANOLEN $NANOPERMPATCH 0.45 1.0 $SPHERE_DIAM $PATCHY_SPHERE_DIAM $ANTIGENE_DIAM 0.0
#create_one_attached <conf_file_name> <nxmax> <nymax> <nzmax> <Lx> <Ly> <Lz> <DensSuperfAntigens> <numspheres> <QFab-diam> <QFab-len> <QFab-diam-permpatch> <QFab-diam-revpatch> <QFab-dist-revpatch> <sphere-diam> <sphere-permpatch-diam> <DiametroAntigene> <bigAntigenSurfDiam> <maxnanobodies>\n
