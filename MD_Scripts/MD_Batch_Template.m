function [TemplateText,qsub,gmx,Server] = MD_Batch_Template(Settings)

if ispc
    Server = 'ced'; % for testing
else
    [~,Servertxt] = system('hostname -s | cut -c 1-3');
    Server = strtrim(Servertxt);
end

if strcmpi(Server,'ced') || strcmpi(Server,'cdr') % cedar
    Account = 'def-thachuk';
elseif ~isempty(regexp(Server,'se[0-9]','ONCE')) || strcmpi(Server,'log') % sockeye
    Account = 'def-thachuk';
else
    Account = 'def-thachuk';
end

if strcmpi(Server,'sea') % Orcinus
    Cores_per_node = 12;
    Nodes = ceil(Settings.Cores/Cores_per_node);
    
    if strcmp(Settings.Mempernode,'-1')
        memline = '';
    elseif strcmp(Settings.Mempernode,'0')
        memline = '';
    else
        memline = ['#PBS -l pmem=' Settings.Mempernode newline];
    end
    
    if Nodes > 1
        mdrun = ['mpirun -np ' num2str(Nodes*2) ...
            ' --npernode $MPI_PPN --mca mpi_paffinity_alone 0 ##MDRUN## -notunepme -ntomp $OMP_NUM_THREADS -maxh ' ...
            num2str(Settings.Hours)];
        gmx = 'gmx_mpi_d';
    else
        mdrun = ['##MDRUN## -maxh ' num2str(Settings.Hours)];
        gmx = 'gmx_d';
    end
    
    TemplateText = ['#!/bin/bash' newline ...
        '#PBS -S /bin/bash' newline ...
        '#PBS -l walltime=' num2str(Settings.Hours) ':' num2str(Settings.Mins) ':00' newline ...
        '#PBS -l nodes=' num2str(Nodes) ':ppn=' num2str(Cores_per_node) newline ...
        memline ...
        '#PBS -V' newline ...
        '#PBS -N ##TASKNAME##' newline ...
        '#PBS -e ##ERROR##.stde' newline ...
        '#PBS -o ##ERROR##.stdo' newline ...
        '#PBS -l partition=QDR' newline ...
        newline newline ...
        '# Check on some basics:' newline ...
        'echo "Running on host: " `hostname`' newline ...
        'echo "Changing to directory from which PBS script was submitted."' newline ...
        'cd ##DIRECTORY##' newline ...
        'echo "Current working directory is now: " `pwd`' newline ...
        newline newline ...
        '# set EXE environment' newline ...
        'module load gromacs/5.1.4' newline ...
        'module load matlab/R2015b' newline ...
        'export MPI_PPN=2' newline ...
        'export OMP_NUM_THREADS=6' newline ...
        newline newline ...
        '# Run Job' newline ...
        '##PREMIN##' newline ...
        mdrun newline ...
        '##CLEANUP##' newline ...
        'echo "Job completed at `date`"' newline ... 
        'exit 0'];

    qsub = 'qsub';
elseif strcmpi(Server,'bel') % Beluga
    Cores_per_node = 40;
    Nodes = ceil(Settings.Cores/Cores_per_node);
    
    if strcmp(Settings.Mempernode,'-1')
        memline = '';
    elseif strcmp(Settings.Mempernode,'0')
        memline = ['#SBATCH --mem-per-cpu=MaxMemPerCPU' newline];
    else
        memline = ['#SBATCH --mem-per-cpu=' Settings.Mempernode newline];
    end
    
    if Nodes > 1
        mdrun = ['srun ##MDRUN## -notunepme -dlb yes -maxh ' num2str(Settings.Hours)];
        gmx = 'gmx_mpi_d';
    else
        mdrun = ['##MDRUN## -maxh ' num2str(Settings.Hours)];
        gmx = 'gmx_d';
    end
    
    TemplateText = ['#!/bin/bash' newline ... 
        '#SBATCH --time=' num2str(Settings.Hours) ':' num2str(Settings.Mins) ':00' newline ... 
        '#SBATCH --nodes=' num2str(Nodes) newline ... 
        '#SBATCH --tasks-per-node=' num2str(Cores_per_node) newline ... 
        memline ...
        '#SBATCH --account=' Account newline ... 
        '#SBATCH --job-name=##TASKNAME##' newline ... 
        '#SBATCH --error=##ERROR##.stde' newline ... 
        '#SBATCH --export=ALL' newline ... 
        newline newline ... 
        '# Check on some basics:' newline ... 
        'echo "Running on host: " `hostname`' newline ... 
        'echo "Changing to directory from which PBS script was submitted."' newline ... 
        'cd ##DIRECTORY##' newline ... 
        'echo "Current working directory is now: " `pwd`' newline ... 
        newline newline ... 
        '# Load modules' newline ... 
        'module load nixpkgs/16.09  gcc/7.3.0  openmpi/3.1.2' newline ... 
        'module load gromacs/2018.3' newline ...
        newline newline ...
        '# Set variables' newline ...
        'export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"' newline ...
        newline newline ...
        '# Run Job' newline ...
        '##PREMIN##' newline ...
        mdrun newline ...
        '##CLEANUP##' newline ...
        'echo "Job completed at `date`"' newline ... 
        'exit 0'];
    
    qsub = 'sbatch';
elseif ~isempty(regexp(Server,'se[0-9]','ONCE')) || strcmpi(Server,'log') % sockeye
    
    Cores_per_node = 32;
    Nodes = ceil(Settings.Cores/Cores_per_node);
    
    if strcmp(Settings.Mempernode,'-1') || strcmp(Settings.Mempernode,'0')
        memline = '186gb';
    else
        memline = Settings.Mempernode;
    end
    
    mdrun = ['mpiexec ##MDRUN## -v -ntomp $OMP_NUM_THREADS -maxh ' ...
        num2str(Settings.Hours)];
    gmx = 'gmx_mpi_d';
    
    TemplateText = ['#!/bin/bash' newline ...
        '#PBS -l walltime=' num2str(Settings.Hours) ':' num2str(Settings.Mins) ':00,' ...
        'select=' num2str(Nodes) ':ncpus=' num2str(Cores_per_node) ':mpiprocs=' num2str(Cores_per_node/8) ...
        ':ompthreads=' num2str(Cores_per_node/4) ':mem=' memline newline ...
        '#PBS -A ' Account newline ...
        '#PBS -N ##TASKNAME##' newline ...
        '#PBS -e ##ERROR##.stde' newline ...
        '#PBS -o ##ERROR##.stdo' newline ...
        newline newline ...
        '# Check on some basics:' newline ...
        'echo "Running on host: " `hostname`' newline ...
        'echo "Changing to directory from which PBS script was submitted."' newline ...
        'cd ##DIRECTORY##' newline ...
        'echo "Current working directory is now: " `pwd`' newline ...
        newline newline ...
        '# set EXE environment' newline ...
        'module load gcc/9.1.0' newline ...
        'module load openmpi/3.1.4' newline ...
        'module load gromacs/5.1.4' newline ...
        'module load matlab/R2018b' newline ...
        newline newline ...
        '# Run Job' newline ...
        '##PREMIN##' newline ...
        mdrun newline ...
        '##CLEANUP##' newline ...
        'echo "Job completed at `date`"' newline ... 
        'exit 0'];

    qsub = 'qsub';
    
elseif strcmpi(Server,'ced') || strcmpi(Server,'cdr') || strcmpi(Server,'gra') % Cedar and graham
    
    if strcmpi(Server,'ced') || strcmpi(Server,'cdr') % Cedar
        Cores_per_node = 48;
    else
        Cores_per_node = 32; % Graham or testing pc
    end
    Nodes = ceil(Settings.Cores/Cores_per_node);
    
    if strcmp(Settings.Mempernode,'-1')
        memline = '';
    elseif strcmp(Settings.Mempernode,'0')
        memline = ['#SBATCH --mem-per-cpu=MaxMemPerCPU' newline];
    else
        memline = ['#SBATCH --mem-per-cpu=' Settings.Mempernode newline];
    end
    
    if Nodes > 1
        mdrun = ['srun ##MDRUN## -notunepme -dlb yes -maxh ' num2str(Settings.Hours)];
        gmx = 'gmx_mpi_d';
    else
        mdrun = ['##MDRUN## -maxh ' num2str(Settings.Hours)];
        gmx = 'gmx_d';
    end
    
    TemplateText = ['#!/bin/bash' newline ... 
        '#SBATCH --time=' num2str(Settings.Hours) ':' num2str(Settings.Mins) ':00' newline ... 
        '#SBATCH --nodes=' num2str(Nodes) newline ... 
        '#SBATCH --tasks-per-node=' num2str(Cores_per_node) newline ... 
        memline ...
        '#SBATCH --account=' Account newline ... 
        '#SBATCH --job-name=##TASKNAME##' newline ... 
        '#SBATCH --output=##ERROR##.stdo' newline ... 
        '#SBATCH --error=##ERROR##.stde' newline ... 
        '#SBATCH --export=ALL' newline ... 
        newline newline ... 
        '# Check on some basics:' newline ... 
        'echo "Running on host: " `hostname`' newline ... 
        'echo "Changing to directory from which PBS script was submitted."' newline ... 
        'cd ##DIRECTORY##' newline ... 
        'echo "Current working directory is now: " `pwd`' newline ... 
        newline newline ... 
        '# Load modules' newline ... 
        'module load nixpkgs/16.09  intel/2016.4 openmpi/2.1.1' newline ... 
        'module load gromacs/2018' newline ...
        'module load matlab/R2018a' newline ...
        newline newline...
        '# Set variables' newline ...
        'export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"' newline ...
        newline newline ...
        '# Run Job' newline ...
        '##PREMIN##' newline ...
        mdrun newline ...
        '##CLEANUP##' newline ...
        'echo "Job completed at `date`"' newline ... 
        'exit 0'];
    
    qsub = 'sbatch';
else
    TemplateText = '';
    qsub = 'local';
    gmx = 'gmx_d';
end

end