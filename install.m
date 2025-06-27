% Installation/Initialization Script for asto-mat
if ~isfolder("./extern/mice")
    rmdir ./extern/mice s
end


%% Specify Versions for Platform
if ispc()
    fprintf("Windows 64 OS Detected\n");
    lsk_string = "latest_leapseconds.tls.pc";

    mice_url = "https://naif.jpl.nasa.gov/pub/naif/toolkit//MATLAB/" + ...
            "PC_Windows_VisualC_MATLAB9.x_64bit/packages/mice.zip";
else
    fprintf("Non-windows architechture\n");
    if ismac()
        fprintf("MAC OS Detected.");
        mice_url = "https://naif.jpl.nasa.gov/pub/naif/toolkit//MATLAB/" + ...
            "MacIntel_OSX_AppleC_MATLAB9.x_64bit/packages/mice.tar.Z";      
    else
        mice_url = "https://naif.jpl.nasa.gov/pub/naif/toolkit//MATLAB" + ...
            "/PC_Linux_GCC_MATLAB9.x_64bit/packages/mice.tar.Z";
    end
    lsk_string = "latest_leapseconds.tls";
end

%% Install MICE
if ispc()
    % Base MICE
    unzip(mice_url, "./extern/");
    else
        command1 = "wget " + mice_url;
        system(command1);
        
        command2 = "tar -xf mice.tar.Z -C extern";
        system(command2);
        
        delete("mice.tar.Z");
    if install_boost
        command1 = "wget " + boost_url;
        system(command1);
    
        command2 = "tar -xf boost_1_87_0.tar.gz";
        system(command2);
    
        delete("boost_1_87_0.tar.gz");
    end
end
fprintf("Successfully installed MICE!\n");

%% Install Kernels

% Save Leapsecond Kernel
websave("data/"+lsk_string, ...
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/" + ...
    lsk_string);
fprintf("Succesfully installed leapsecond kernel!\n");

% Save Planetary Ephemerides
websave("data/de440.bsp", "https://naif.jpl.nasa.gov/pub/naif/" + ...
    "generic_kernels/spk/planets/de440.bsp");
fprintf("Successfully added DE440\n");

% Save GM Data
websave("data/gm_de440.tpc", "https://naif.jpl.nasa.gov/pub" + ...
    "/naif/generic_kernels/pck/gm_de440.tpc");

%% Path Configuration
addpath(genpath("./"));

savepath([userpath,filesep,'pathdef.m']);