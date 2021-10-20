Hi Whale scientists!

Those tools have been developed to facilitate and accelerate the beamforming algorithm for the BRING project.

In the model, you will get the information of the array automatically by specifying the arrID which can take 4 different values: AAV/CLD/MLB/PRC

The best is to start from the config folder and do a copy of the defautlConfig.m script and rename it for your actual config run. I suggest starting a new folder (ex. configKevMac) for each user so we do not have any trouble with path folder.

Essentially the information needed to run the script are in the defautlConfig.m file. In this file, you will set all the parameters needed for the script to operate. Then you can add figures, more data treatment and all kind of things after the call of locateBring.m . The main idea is to have the configXX.m file to set all information about the data to open and the parameter to use, keep the locateBring.m just a hard code machine that does the math, some showXXX.m script that plot figure and do some post-treatment.

There is some script variation. For quick use without opening and copying other script you can use directly the locateScript.m but this might not be updated. Another file that can be used is the locateFun. This function is call as follow:

[angleM, matEnergie] = locateFun(arrID, ptime, varargin).

Where the varargin option available like Ns, Fs_min, Fs_max can be specified and are list in the locateInput.m file. otherwise, default parameter will be used. This can be useful for a quick plot, treatment or to get the value of a sound in a forloop.

Your folder three should be {Yourpath}/{arrID}/ so each folder will contain the WAV related to one array. Don't forget you «/» after folder name as the code will add outName to it to complete the total path.

You will need m_map and probably better to get the gebco_colormap.dat and the GSL_bathy_500m.mat file and add this to your permanent matlab path addpath({path_to_map_ool}). Those file need will probably be send to you by email. Otherwise email me or look online.

-----------------------------------
For your first use of the tool please follow these step:
1) Clone rep from git hub : git@github.com:duik15/Localisation.git
2) Install m_map, better get the gebco_colormap.dat and the GSL_bathy_500m.mat file
and add this to your permanent matlab path addpath({path_to_map_ool})
3) Create a config folder( mkdir configKev) inside localisation folder
4) Copy and paste the defaultConfig.m in oyur new configFolder
5) Modify the output and input path in the scrip
6) Create an output tree folder as follow {yourPreferedPath}/results/arrID{AAV/CLD...}
7) Modifiy location/.gitignore to add your configFolder to ignore git tracking
8) Open localisation/getTools/getMachine.m and add your computer name so you can use this tool. Now you can get specify path related to your computer using various other get function like loadGPXTrack.m.

Enjoy!
----------------------------------

Do not hesitate to participate in the development of those tools. You can email me at kevin.duquette@dfo-mpo.gc.ca for any questions.
