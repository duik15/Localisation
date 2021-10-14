Hi Whale scientists! 

Those tools have been developped to facilite and accelerate the beamforming algorithme for the BRING projet. 

The model you will the information of the array automatically by specified the arrID which can take 4 different value: AAV/CLD/MLB/PRC

The best is to start from the config folder and do a copy of the defautlConfig.m script and rename it for your actual config run. I suggest starting a new folder (ex. configKevMac) for each user so we do not have any trouble with path folder. 

Ensentially the information needed to run the script are in the defautlConfig.m file. In this file you will set all the parameters needed for the script to operate. Then you can add figures, more data treatment and all kind of things after the call of locateBring.m . The main idea is to have the configXX.m file to set all information about the data to open and the parameter to use, keep the locateBring.m just a hard code machine that do the math, some showXXX.m script that plot figure and do some post treatment. 

There is some script variation. For a quick use without opening and copying other scritpt you can use directly the locateScript.m but this might not be updated. Another file that can be use is the locateFun. This function is call as follow:

[angleM, matEnergie] = locateFun(arrID, ptime, varargin).

Where the varargin otpion available like Ns, Fs_min, Fs_max can be specified and are list in the locateInput.m file. otherwise, default parameter will be use. This can be usefull for a quick plot, treatment or to get the value of a sound in a forloop.

Your folder three should be {Yourpath}/{arrID}/ so each folder will contain the WAV related to one array. Don't forget you «/» after folder name as the code will add outName to it to complet the total path. 

Do not hesitate to participate in the developpement of those tools. You can email me at kevin.duquette@dfo-mpo.gc.ca for any questions. 
