function [ phiHO, ephiHO, jHO ] = extractSource(phi_hdf5file, j_hdf5file, data)


hinfo = hdf5info(phi_hdf5file);

phiHO = hdf5read(hinfo.GroupHierarchy.Datasets(2));

%phiHO = flipud(phiHO); % UNSORTED VERSION OF PROFUGUS, DELETE WHEN UPDATED.


hinfo = hdf5info(j_hdf5file);

jHO = hdf5read(hinfo.GroupHierarchy.Datasets(2));


ephiHO = hdf5read(hinfo.GroupHierarchy.Datasets(5)); % 5th item in the inf_med_current.h5 set


end

