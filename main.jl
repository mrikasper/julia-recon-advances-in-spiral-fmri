#############################################
## Set parameters
#############################################

using NIfTI, PyPlot, HDF5, MRIReco, LinearAlgebra, Dates

# parameters
@info "Setting Parameters"

do_recalc_sensitivity = false
do_b0_correction = true
do_inspect_iterations = false

# from previous Matlab calculation
# rad/s, manually determined from conversion of b0 map
w0_offset = + 76.4903
dt = 1.8e-6 # acquisition dwell time [s]
TE = 0.02 # s
Nx = 300 # 240
Ny = 365 # 292
fov = [190, 230, 0.9+0.1]/1000 # [m]

idx_slice = 18; # slice(s) to reconstruct
do_select_slice =  !isempty(idx_slice)
# Paths
path_data = "D:\\SPIFI\\ExportETHResearchCollection"

# Raw Coil/Trajectory data: ISMRMRD file (.nii); units: rad/m
filename_data = joinpath(path_data,"SPIFI_0007_RawData_RotatedToSliceGeometry3D_spiralOut_singleVolume.h5")

# Maps are NIFTI files (.nii); units: Hz
filename_sense_magnitude = joinpath(path_data, "SPIFI_0007_MapsForReconstruction", "coilSensitivityMaps_magnitude.nii")
filename_sense_phase = joinpath(path_data, "SPIFI_0007_MapsForReconstruction", "coilSensitivityMaps_phase.nii")
filename_b0 = joinpath(path_data, "SPIFI_0007_MapsForReconstruction", "b0MapSmoothedResizedSpiralGeometry1_Hz.nii")

# rawdata: ch, il, read was order
# traj: il, read, kdim(x,y,z)
# Julia order: read, il, ch and kdim, read, il
#data = permutedims(h5read(filename_data, "rawdata"),[3,2,1,4])


#############################################
## Load ISMRMRD (coil/traj) data and maps (SENSE/B0)
#############################################

@info "Load Data"
dataFile = ISMRMRDFile(filename_data)
rawData = RawAcquisitionData(dataFile)

# select data from selected slices only
if do_select_slice
	indices = 1:length(rawData.profiles)
	#ic = [x for x ∈ indices if x ∉ idx_slice]
	ic = setdiff(indices, idx_slice)
	deleteat!(rawData.profiles,ic)

	n_slices =  size(idx_slice,1);
	# make sinlge slice parameter set
	rawData.profiles[1].head.idx.slice = 0;
	rawData.params["enc_lim_slice"] = Limit(0, 0, 0)
	rawData.params["reconSize"][3] = n_slices
	rawData.params["encodedSize"][3] = n_slices
	rawData.params["reconFOV"][3] = rawData.params["reconFOV"][3]/size(indices,1)*n_slices
	rawData.params["encodedFOV"][3] = rawData.params["encodedFOV"][3]/size(indices,1)*n_slices
	# print function for all parameters
	 for (key, value) in rawData.params; print(key); print(": ");print(value);print("\n"); end
end

@info "Converting RawAcquisitionData to AcquisitionData"
# acqData = AcquisitionData(rawData,estimateProfileCenter=false)
acqData = AcquisitionData(rawData)

# adjust TE and TAQ after read-in, is not taken from file
n_samples = rawData.params["encodedSize"][1]
tAQ = (n_samples-1) * dt
acqData.traj[1].AQ=tAQ # important for B0 correction
acqData.traj[1].TE= TE
acqData.traj[1].times = TE .+ collect(0:dt:tAQ)
#############################################
## Load and convert traj from rad/m to -0.5 -> 0.5
#############################################

@info "Normalize trajectory rad/m -> [-1/2 1/2] FOV"
traj_node_normalized = acqData.traj[1].nodes
traj_node_normalized /= (2*π)
traj_node_normalized[2,:,:] /= Nx/fov[1]
traj_node_normalized[1,:,:] /= Ny/fov[2]

# TODO: remove 3rd dimension for recon, otherwise MRIReco.jl has 3D assumption
traj_node_normalized = traj_node_normalized[1:2,:]

acqData.traj[1].nodes = traj_node_normalized


#############################################
## Plot Trajectory
##############################################

@info "Plot Normalized Trajectory"

mytr = acqData.traj[1]
figure(1);
subplot(1,2,2)
plot(mytr.nodes[1,:], mytr.nodes[2,:]);
#plot(traj[1,:], traj[2,:]);
axis(:square)
gcf()

################################
## Generate coil sensitivity maps
################################

if do_recalc_sensitivity == true
    @info "Calc Coil Sensitivity maps with Espirit"
    acqDataCart = regrid2d(acqData, (Nx,Ny); cgnr_iter=3)
    sense = espirit(acqDataCart,(6,6),30,eigThresh_1=0.02, eigThresh_2=0.98)
    sensitivity = sense;
#    sensitivity[1,1,1] = reshape(sense,:,n_channels)
else
    @info "Load Coil Sensitivity maps"
    # from previous Matlab calculation

	#loadImage(filename_sense_magnitude)
	sense = niread(filename_sense_magnitude).*exp.(1im.*niread(filename_sense_phase))

	n_channels = size(sense,4);

	sense = sense[:,:,idx_slice,:]
	# keep channels in 4th dim, even if singleton slice dim
	if ndims(sense)==3
		size_sense = size(sense)
		sense = reshape(sense, size_sense[1], size_sense[2], 1, size_sense[3])
	end

	n_slices = size(sense,3);

    sensitivity = Array{Array{Complex{Float64},2},4}(undef,1,1,1,1)
    sensitivity = sense;
end

sensitivity = mapslices(x ->imresize(x, (Nx, Ny)), sensitivity, dims=[1,2])
# conversion to ComplexF64 needed for reconstruction_2d solver
sensitivity = convert(Array{ComplexF64, 4}, reshape(sensitivity, Nx, Ny, 1, n_channels))


figure(2); clf(); for ch in 1:n_channels; subplot(8,4,ch); imshow(rotl90(abs.(sense[:,:,1,ch]))); end;
subplots_adjust(wspace=0.05,hspace=0.05,left=0.05,bottom=0.0,right=1.0,top=0.95)
gcf()
figure(3); clf(); for ch in 1:n_channels; subplot(8,4,ch); imshow(rotl90(angle.(sense[:,:,1,ch]))); end;
subplots_adjust(wspace=0.05,hspace=0.05,left=0.05,bottom=0.0,right=1.0,top=0.95)
gcf()

##########################
## Load B0 map, adapt geometry
##########################

@info "Load B0 map"
b0 = niread(filename_b0).*2π .- w0_offset

#select recon slice
b0 = b0[:,:,idx_slice]

resizedB0 = mapslices(x->imresize(x,(Nx, Ny)), b0, dims=[1,2])

cmap = 1im.*resizedB0;

figure(4); cla(); imshow(rotl90(resizedB0), cmap="gray");
subplots_adjust(wspace=0.05,hspace=0.05,left=0.05,bottom=0.0,right=1.0,top=0.95)
colorbar()
gcf()


##########################
## Perform reference reconstruction
##########################
@info "Reference cg-SENSE recon"
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (Nx,Ny)
params[:regularization] = "L2"
params[:λ] = 1.e-2
params[:iterations] = 10
params[:solver] = "cgnr"
params[:solverInfo] = SolverInfo(ComplexF64,store_solutions=do_inspect_iterations)
params[:senseMaps] = sensitivity

# params[:reco] = "direct"
@time begin
    if do_b0_correction
        params[:correctionMap] = cmap
        #params[:alpha] = 1.75 # oversampling factor for interpolation
        #params[:m] = 4.0 # truncation size of interpolation kernel
        #params[:K] = 28  # number of translates for LeastSquares approaches (not NFFT-approach) to time-segmentation
    end

    img_ref = reconstruction(acqData, params).data
end

# h5write(filename_recon, "/img_ref", img_ref)


##########################
## Plot resulting recon and iterations
##########################

figure(5); cla();
subplot(1,2,1); imshow(rotl90(abs.(img_ref[:,:,1,1,1])), cmap="gray", vmax=0.8, aspect="equal");
subplot(1,2,2); imshow(rotl90(angle.(img_ref[:,:,1,1,1])), cmap="gray");
subplots_adjust(wspace=0.05,hspace=0.05,left=0.05,bottom=0.0,right=1.0,top=0.95)
gcf()

figure(6); cla();
subplot(1,2,1); plot(params[:solverInfo].convMeas); title("Linear")
subplot(1,2,2); semilogy(params[:solverInfo].convMeas);  title("Log")
gcf().suptitle("Convergence over Iterations")
xlabel("Iterations")
gcf()


# plot with colorbar
# begin
# 	fig, (ax1, ax2) = subplots(figsize=(9, 3), ncols=2)
# 	hp1 = ax1.imshow(rotl90(abs.(img_ref[:,:,1,1,1])), cmap="gray", vmax=0.8, aspect="equal")
# 	hp2 = ax2.imshow(rotl90(angle.(img_ref[:,:,1,1,1].-π/10.0)), cmap="gray")
# 	fig.colorbar(hp1, ax=ax1)
# 	fig.colorbar(hp2, ax=ax2)
# 	fig
# end

# Plot Iterations
if do_inspect_iterations
 	img_iter = Vector{Array{ComplexF64,5}}(undef,11)
 	for iter = 1:11
     	img_iter[iter] = reshape(params[:solverInfo].x_iter[iter],Nx,Ny,1,1,1)
     	figure(1000+iter);cla;imshow(abs.(img_ref[:,:,1,1,1])); gcf()
	end
end

##########################
## Save final recon result with a time stamp
##########################
figure(5);gcf()
figSaveName = string("FinalRecon_", Dates.format(Dates.now(), "yyyy-mm-dd_HH_MM_SS"), ".png")
savefig(figSaveName,dpi=300)
