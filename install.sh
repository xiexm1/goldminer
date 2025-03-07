#!/usr/bin/env sh
stat=0

echo "check software dependencies"

if [ `which mcl | wc -l` == 0 ];then
	echo 'ERROR: mcl cannot be found'
	stat=$((stat+1))
fi

if [ `which R | wc -l` == 0 ];then
	echo 'ERROR: R cannot be found'
	stat=$((stat+1))
fi

if [ `which python | wc -l` == 0 ];then
	echo 'ERROR: python cannot be found'
	stat=$((stat+1))
fi

if [ $stat != 0 ];then
	exit 1
fi

rm -rf bin/
mkdir bin/
cd bin/
	for ph in ../src/*R ../src/*py ../src/*sh;do
		pc=`basename $ph | cut -d"." -f1`
		echo "bulding "`basename $ph`
		ln -s $ph $pc
		chmod +x $pc
	done
cd ..

echo ""
echo "Installation in finished"
echo ""
echo "Please add the following line to your ~/.bash_profile, and source ~/.bash_profile before running goldminer."
echo ""
echo '	export PATH='`pwd`':$PATH'
echo ""
