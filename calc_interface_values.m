function [psix, psiy] = calc_interface_values(psi)
%  Calculate interface values for the matrix interfaces
%  ---------------------------------------------------------------------
%  Copyright (C) 2016 by the LearnEDFM authors
% 
%  This file is part of LearnEDFM.
% 
%  LearnEDFM is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  LearnEDFM is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with LearnEDFM.  If not, see <http://www.gnu.org/licenses/>.
%  ---------------------------------------------------------------------
% 
%  Authors: Gunnar Jansen, University of Neuchatel, 2016 
%
%  calc_frac_mean_mat(N_fractures,Nf_f, Nf_i)
%
%  Input: 

%  Output:
% 

Nf = size(psi);

psix             = sparse(Nf(1)+1,Nf(2));                                                                                                                
psiy             = sparse(Nf(1),Nf(2)+1);
psix(2:Nf(1),:)  = 2./(1./psi(1:Nf(1)-1,:) + 1./psi(2:Nf(1),:));             
psiy(:,2:Nf(2))  = 2./(1./psi(:,1:Nf(2)-1) + 1./psi(:,2:Nf(2)));
psix(1,:)        = psi(1,:);                                    
psix(Nf(1)+1,:)  = psi(Nf(1),:);
psiy(:,1)        = psi(:,1);
psiy(:,Nf(2)+1)  = psi(:,Nf(2));