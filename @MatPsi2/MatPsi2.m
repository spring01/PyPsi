classdef MatPsi2 < handle
    
    properties (SetAccess = private)
    
        pathMatPsi2; % Path of @MatPsi2 folder
        
    end
    
    properties (Access = private, Transient = true)
    
        objectHandle; % Handle to the underlying C++ class instance
        
    end
    
    methods
        %% Constructor - Create a new C++ class instance  
        function this = MatPsi2(molStr, basisSet, charge, multiplicity, psiDataDir)
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
            this.objectHandle = MatPsi2.MatPsi2_mex('new', molStr, basisSet, charge, multiplicity, psiDataDir);
        end
        
        %% Destructor - Destroy the C++ class instance 
        function delete(this)
            if(~isempty(this.objectHandle))
                MatPsi2.MatPsi2_mex('delete', this.objectHandle);
            end
        end
        
        function varargout = InputInfo_MoleculeString(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('InputInfo_MoleculeString', this.objectHandle, varargin{:});
        end
        
        function varargout = InputInfo_BasisSet(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('InputInfo_BasisSet', this.objectHandle, varargin{:});
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
        
        function varargout = Molecule_NuclearRepulsionEnergy(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Molecule_NuclearRepulsionEnergy', this.objectHandle, varargin{:});
        end
        function varargout = Molecule_SetCharge(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Molecule_SetCharge', this.objectHandle, varargin{:});
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
        
        function varargout = BasisSet_FunctionToCenter(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_FunctionToCenter', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_FunctionToAngularMomentum(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_FunctionToAngularMomentum', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_PrimitiveExponents(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_PrimitiveExponents', this.objectHandle, varargin{:});
        end
        
        function varargout = BasisSet_PrimitiveCoefficients(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('BasisSet_PrimitiveCoefficients', this.objectHandle, varargin{:});
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
        
        function varargout = Integrals_PotentialPointCharges(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Integrals_PotentialPointCharges', this.objectHandle, varargin{:});
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
        
        function varargout = Integrals_IndicesForExchange(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('Integrals_IndicesForExchange', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_Initialize(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_Initialize', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_Type(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_Type', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_DensityToJ(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_DensityToJ', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_DensityToK(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_DensityToK', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_OrbitalToJ(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_OrbitalToJ', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_OrbitalToK(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_OrbitalToK', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_OccupiedOrbitalToJ(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_OccupiedOrbitalToJ', this.objectHandle, varargin{:});
        end
        
        function varargout = JK_OccupiedOrbitalToK(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('JK_OccupiedOrbitalToK', this.objectHandle, varargin{:});
        end
        
        function varargout = DFJK_mnQMatrixUnique(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('DFJK_mnQMatrixUnique', this.objectHandle, varargin{:});
        end
        
        function varargout = DFJK_mnQTensorFull(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('DFJK_mnQTensorFull', this.objectHandle, varargin{:});
        end
        
        function varargout = DFJK_mnAMatrixUnique(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('DFJK_mnAMatrixUnique', this.objectHandle, varargin{:});
        end
        
        function varargout = DFJK_mnATensorFull(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('DFJK_mnATensorFull', this.objectHandle, varargin{:});
        end
        
        function varargout = DFJK_InverseJHalfMetric(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('DFJK_InverseJHalfMetric', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_DoSCF(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_DoSCF', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_Reset(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_Reset', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_EnableMOM(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_EnableMOM', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_DisableMOM(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_DisableMOM', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_EnableDamping(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_EnableDamping', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_isableDamping(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_DisableDamping', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_EnableDIIS(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_EnableDIIS', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_DisableDIIS(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_DisableDIIS', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_GuessSAD(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_GuessSAD', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_GuessCore(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_GuessCore', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_TotalEnergy(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_TotalEnergy', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_Orbital(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_Orbital', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_OrbitalEnergies(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_OrbitalEnergies', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_Density(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_Density', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_CoreHamiltonian(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_CoreHamiltonian', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_JMatrix(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_JMatrix', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_KMatrix(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_KMatrix', this.objectHandle, varargin{:});
        end
        
        function varargout = RHF_FockMatrix(this, varargin)
            [varargout{1:nargout}] = MatPsi2.MatPsi2_mex('RHF_FockMatrix', this.objectHandle, varargin{:});
        end

    end
    
    methods (Static, Access = private)
        
        [varargout] = MatPsi2_mex(command_name, objectHandle, varargin);
        
    end
    
end
