all:
	cd my_routines; make
	cd arpack; make
	cd working; make
	
unbackup:
	untar -xvf all.tar.gz
	cd my_routines; tar -xvf my_routines.tar.gz
	cd arpack; tar -xvf arpack.tar.gz
	cd working; tar -xvf working.tar.gz

backup:
	cd my_routines; make backup
	cd arpack; make backup
	cd working; make backup
	tar --ignore-failed-read -cvzf others.tar.gz ./working/data/*.cpp ./working/data/*.h 
	tar --ignore-failed-read -cvzf all.tar.gz ./my_routines/my_routines.tar.gz ./arpack/arpack.tar.gz ./working/working.tar.gz others.tar.gz makefiles makefile

clean:
	cd my_routines; make clean
	cd arpack; make clean
	cd working; make clean

server-upload:
	make backup
	scp all.tar.gz pstano@quantum.quniverse.sk:
	ssh pstano@quantum.quniverse.sk "./quniverse2server $1 $2"

server-download:
	ssh pstano@quantum.quniverse.sk "./server2quniverse $1 $2"
	scp pstano@quantum.quniverse.sk:results.tar.gz .
	
quniverse2server:
	scp -o BindAddress=147.213.112.193 all.tar.gz pstano@192.168.1.20:CPP
	echo $1 will be run with the result directory data/$2
	ssh -b 147.213.112.193 pstano@192.168.1.20 "cd CPP; make unbackup; rm -f makefile_settings; ln -s makefiles/makefile_settings_virtual_server makefile_settings; make clean; make; cd working; mkdir data/$2; nohup $1 > data/$2/terminal_output 2> data/$2/terminal_output"  
	
server2quniverse:
	ssh -b 147.213.112.193 pstano@192.168.1.20 "tar -cvf results.tar.gz CPP/working/data/$2"
	scp -o BindAddress=147.213.112.193 pstano@192.168.1.20:results.tar.gz .
	
	
	


