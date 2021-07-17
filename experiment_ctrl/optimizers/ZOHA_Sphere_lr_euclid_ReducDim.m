classdef ZOHA_Sphere_lr_euclid_ReducDim < ZOHA_Sphere_lr_euclid
    properties
        code_len % dimension of the full space
        % N % dimension of the subspace
        basis
    end
    
    methods
        function self = ZOHA_Sphere_lr_euclid_ReducDim(space_dimen, subspac_d, options)
            self@ZOHA_Sphere_lr_euclid(subspac_d, options);
            self.code_len = space_dimen;
            % basis = self.getBasis("rand");
            self.opts.subspac_d = self.dimen;
        end
        
        function basis = getBasis(self, basis_opt, varargin)
            % basis could be a vector set or a path to some mat file or
            % random direction
            if isstring(basis_opt) || ischar(basis_opt)
                if basis_opt == "rand" % use rand N dimension (rnd orthogonal basis) in the space 
                    subspacedim = self.dimen; %varargin(1); %  
                    code_length = self.code_len;
                    basis = zeros(subspacedim, code_length);
                    for i = 1:subspacedim
                        tmp_code = randn(1, code_length);
                        tmp_code = tmp_code - (tmp_code * basis') * basis;
                        basis(i, :) = tmp_code / norm(tmp_code);
                    end
                else
                    try % load the file from 
                        basis = load(basis_opt, 'PC_codes');
                    catch ME
                        error("%s happend when trying to load the PC_codes.",ME);
                    end
                end
            else
                assert(isa(basis_opt,'double'), 'Unrecoginized data type for basis option.')
                basis = basis_opt;
            end
            if all(size(basis) == [self.dimen, self.code_len])
            elseif all(size(basis) == [self.code_len, self.dimen])
                basis = basis';
            else
                error("Size of imported basis (%d,%d) doesn't match that predefined in the Optimizer (%d,%d).", ...
                    size(basis,1), size(basis,2), self.dimen, self.code_len)
            end
            self.basis = basis;
        end
        
        function parseParameters(self, opts)
            self.parseParameters@ZOHA_Sphere_lr_euclid(opts);
        end
        function [new_samples, new_ids] =  doScoring(self, codes, scores, maximize, TrialRecord)
            codes_coord = codes * self.basis'; % isometric transform by the basis set
            [new_coord, new_ids] = self.doScoring@ZOHA_Sphere_lr_euclid(codes_coord, scores, maximize, TrialRecord);
            new_samples = new_coord * self.basis; % isometric inverse transform by the basis set
        end
    end
end