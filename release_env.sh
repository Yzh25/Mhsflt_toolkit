scriptdir=`pwd`
zumisenv=${scriptdir}/zUMIs-env
miniconda=${scriptdir}/zUMIs-miniconda.tar.bz2

#relase report module
cat ${scriptdir}/report/create_summary.parta* > ${scriptdir}/report/create_summary
chmod 755 ${scriptdir}/report/create_summary

#check if zUMIs environment has been unpacked from tar
if [[ ! -d ${zumisenv} ]] || [[ ${zumisdir}/zUMIs-miniconda.partaa -nt ${zumisenv} ]] ; then
    [ -d ${zumisenv} ] || mkdir -p ${zumisenv}
    cat `pwd`/zUMIs-miniconda.parta* > ${miniconda}
    tar -xj --overwrite -f ${miniconda} -C ${zumisenv}
fi
#activate zUMIs environment!
  
unset PYTHONPATH
unset PYTHONHOME
source ${zumisenv}/bin/activate
conda-unpack

echo "environment release finish"
