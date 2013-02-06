while true; do
    read -p "Do you wish to initialize sus simulation killing all jobs? " yn
    case $yn in
        [Yy]* ) echo "ok, initializing..."; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done 
if [ -e MOS_JOB_ID ]
then
moskillall  -J`cat MOS_JOB_ID`
rm MOS_JOB_ID
fi
exit
while true; do
    read -p "Do you wish to delete all directories? " yn
    case $yn in
        [Yy]* ) echo "ok, deleting..."; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done 
rm -fr N_*_*
