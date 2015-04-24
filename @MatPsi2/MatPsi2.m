classdef MatPsi2 < handle
    
    properties (SetAccess = private)
    
        pathMatPsi2; % Path of @MatPsi2 folder
        
    end
    
    properties (Access = private, Transient = true)
    
        objectHandle; % Handle to the underlying C++ class instance
        
    end
    
    methods
        %% Constructor - Create a new C++ class instance  
        function this = MatPsi2(cartesian, basisSet, charge, multiplicity, psiDataDir)
            if(nargin < 3)
                charge = 0;
            end
            if(nargin < 4)
                multiplicity = 1;
            end
            if(exist('./@MatPsi2', 'file'))
                pathMatPsi2 = [pwd(), '/@MatPsi2'] ;
            else
                currpath = path();
                paths_num = length(regexp(currpath, ':', 'match')) + 1;
                for i = 1:paths_num
                    toppath = regexp(currpath, '[^:]*', 'match', 'once');
                    if(exist([toppath, '/@MatPsi2'], 'file'))
                        pathMatPsi2 = [toppath, '/@MatPsi2'];
                        break;
                    else
                        currpath = currpath(length(toppath)+2:end);
                    end
                end
                if(i>=paths_num)
                    disp('MatPsi2 cannot find itself; an exception might be thrown soon.');
                    pathMatPsi2 = [];
                end
            end
            if(nargin < 5)
                psiDataDir = pathMatPsi2;
            end
            this.pathMatPsi2 = pathMatPsi2;
            this.objectHandle = MatPsi2.MatPsi2_mex('new', cartesian, basisSet, charge, multiplicity, psiDataDir);
        end
        
        %% Destructor - Destroy the C++ class instance 
        function delete(this)
            if(~isempty(this.objectHandle))
                MatPsi2.MatPsi2_mex('delete', this.objectHandle);
            end
        end
        
        function varargout = Settings_MaxNumCPUCores(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Settings_MaxNumCPUCores', this.objectHandle, varargin{:});
        end
        
        function varargout = Settings_MaxMemoryInGB(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Settings_MaxMemoryInGB', this.objectHandle, varargin{:});
        end
        
        function varargout = Settings_PsiDataDir(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Settings_PsiDataDir', this.objectHandle, varargin{:});
        end
        
        function varargout = Settings_TempDir(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Settings_TempDir', this.objectHandle, varargin{:});
        end
        
        function varargout = Settings_SetMaxNumCPUCores(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Settings_SetMaxNumCPUCores', this.objectHandle, varargin{:});
        end
        
        function varargout = Settings_SetMaxMemory(this, varargin)
            if(isfloat(varargin{1}))
                varargin{1} = num2str(varargin{1});
            end
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Settings_SetMaxMemory', this.objectHandle, varargin{:});
        end
        
        function varargout = Settings_SetPsiDataDir(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Settings_SetPsiDataDir', this.objectHandle, varargin{:});
        end
        
        function varargout = Molecule_Fix(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Molecule_Fix', this.objectHandle, varargin{:});
        end
        
        function varargout = Molecule_Free(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Molecule_Free', this.objectHandle, varargin{:});
        end
        
        function varargout = Molecule_NumAtoms(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Molecule_NumAtoms', this.objectHandle, varargin{:});
        end
        
        function varargout = Molecule_NumElectrons(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Molecule_NumElectrons', this.objectHandle, varargin{:});
        end
        
        function varargout = Molecule_Geometry(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Molecule_Geometry', this.objectHandle, varargin{:});
        end
        
        function varargout = Molecule_SetGeometry(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Molecule_SetGeometry', this.objectHandle, varargin{:});
        end
        
        function varargout = Molecule_AtomicNumbers(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Molecule_AtomicNumbers', this.objectHandle, varargin{:});
        end
        
        function varargout = Molecule_NucRepEnergy(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Molecule_NucRepEnergy', this.objectHandle, varargin{:});
        end
        function varargout = Molecule_SetCharge(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Molecule_SetCharge', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_Name(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_Name', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_SetBasisSet(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_SetBasisSet', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_IsSpherical(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_IsSpherical', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_NumFunctions(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_NumFunctions', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_NumShells(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_NumShells', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_ShellTypes(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_ShellTypes', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_ShellNumPrimitives(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_ShellNumPrimitives', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_ShellNumFunctions(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_ShellNumFunctions', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_ShellToCenter(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_ShellToCenter', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_FuncToCenter(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_FuncToCenter', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_FuncToShell(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_FuncToShell', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_FuncToAngular(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_FuncToAngular', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_PrimExp(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_PrimExp', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_PrimCoeffUnnorm(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_PrimCoeffUnnorm', this.objectHandle, varargin{:});
        end
        
        function varargout = Integrals_Overlap(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Integrals_Overlap', this.objectHandle, varargin{:});
        end
        
        function varargout = Integrals_Kinetic(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Integrals_Kinetic', this.objectHandle, varargin{:});
        end
        
        function varargout = Integrals_Potential(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Integrals_Potential', this.objectHandle, varargin{:});
        end
        
        function varargout = Integrals_PotentialEachCore(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Integrals_PotentialEachCore', this.objectHandle, varargin{:});
        end
        
        function varargout = Integrals_PotentialPtQ(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Integrals_PotentialPtQ', this.objectHandle, varargin{:});
        end
        
        function varargout = Integrals_Dipole(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Integrals_Dipole', this.objectHandle, varargin{:});
        end
        
        function varargout = Integrals_ijkl(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Integrals_ijkl', this.objectHandle, varargin{:});
        end
        
        function varargout = Integrals_NumUniqueTEIs(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Integrals_NumUniqueTEIs', this.objectHandle, varargin{:});
        end
        
        function varargout = Integrals_AllUniqueTEIs(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Integrals_AllUniqueTEIs', this.objectHandle, varargin{:});
        end
        
        function varargout = Integrals_AllTEIs(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Integrals_AllTEIs', this.objectHandle, varargin{:});
        end
        
        function varargout = Integrals_IndicesForK(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Integrals_IndicesForK', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_Initialize(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_Initialize', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_Type(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_Type', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_DensToJ(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_DensToJ', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_DensToK(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_DensToK', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_OrbToJ(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_OrbToJ', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_OrbToK(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_OrbToK', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_OccOrbToJ(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_OccOrbToJ', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_OccOrbToK(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_OccOrbToK', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_DFTensor_AuxPriPairs(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_DFTensor_AuxPriPairs', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_DFTensor_AuxPriPri(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_DFTensor_AuxPriPri', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_DFMetric_InvJHalf(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_DFMetric_InvJHalf', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_RunRHF(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_RunRHF', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_RunUHF(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_RunUHF', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_RunRKS(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_RunRKS', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_EnableMOM(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_EnableMOM', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_DisableMOM(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_DisableMOM', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_EnableDamping(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_EnableDamping', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_DisableDamping(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_DisableDamping', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_EnableDIIS(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_EnableDIIS', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_DisableDIIS(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_DisableDIIS', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_GuessSAD(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_GuessSAD', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_GuessCore(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_GuessCore', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_TotalEnergy(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_TotalEnergy', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_OrbitalAlpha(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_OrbitalAlpha', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_OrbitalBeta(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_OrbitalBeta', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_OrbEigValAlpha(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_OrbEigValAlpha', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_OrbEigValBeta(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_OrbEigValBeta', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_DensityAlpha(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_DensityAlpha', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_DensityBeta(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_DensityBeta', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_CoreHamiltonian(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_CoreHamiltonian', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_FockAlpha(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_FockAlpha', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_FockBeta(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_FockBeta', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_Gradient(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_Gradient', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_InitialGuessDensity(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_InitialGuessDensity', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_RHF_Coulomb(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_RHF_Coulomb', this.objectHandle, varargin{:});
        end
        
        function varargout = SCF_RHF_Exchange(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('SCF_RHF_Exchange', this.objectHandle, varargin{:});
        end

    end
    
    methods (Static, Access = private)
        
        [varargout] = MatPsi2_mex(command_name, objectHandle, varargin);
        
    end
    
end
