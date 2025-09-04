function [Poly,InvCheb,Gamma]=createSparsebasis(grid,maxdim,Xss)
% This Function creates a sparse Basis, POLY, for polynomial approximation.
% It assumes that the grid points coincide with the Chebyshev nodes for a
% relevant transformation of the state space/ basis functions.
% It also provides the generalized inverse, INVCHEB, to estimate
% coefficients of the polynomial from function valuies using least squares.
%
% In addition it provides a matrix GAMMA such that GAMMA*x is a
% perturbation of the marginal distribution functions stored in XSS
%
% Authors: Christian Bayer and Ralph C. Luetticke, Bonn and London.
%
% =========================================================================
% Part of the Matlab code to accompany the paper
% 'Precautionary Savings, Illiquid Assets, and the Aggregate Consequences of
%  Shocks to Household Income Risk', Bonn mimeo
% http://wiwi.uni-bonn.de/hump/wp.html
% =========================================================================
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The most recent version or successor of the above-mentioned paper is
% properly cited in all work and publications that benefit from insights
% derived from the Software.
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%
% =========================================================================


Poly=[];
maxlevel=max([grid.nb grid.na grid.nse]);

%% Values of Chebyshev Polynomials (1st type) at Chebyshev Nodes

Tb  = cos(pi* (0:maxlevel-1)' * (linspace(0.5/grid.nb/2,1-0.5/grid.nb*2,grid.nb)))';
Ta  = cos(pi* (0:maxlevel-1)' * (linspace(0.5/grid.na*2,1-0.5/grid.na/2,grid.na)))';
Tse = cos(pi* (0:maxlevel-1)' * (linspace(0.5/(grid.nse-1),1-0.5/(grid.nse-1),(grid.nse-1))))';


for j1=1:(length(grid.se)-1)
    for j2=1:length(grid.a)
        for j3=1:length(grid.b)
            if j1+j2+j3<maxdim
                [TT1, TT2, TT3] = ndgrid(Tb(:,j3), Ta(:,j2), [Tse(:,j1); 0]);
                Poly=[Poly TT1(:).*TT2(:).*TT3(:)];
            end
        end
    end
end

for j2=1:length(grid.b)
    for j3=1:length(grid.a)
        if j2+j3<maxdim-1
            [TT1, TT2, TT3] = ndgrid(Tb(:,j2), Ta(:,j3), [zeros(length(grid.se)-1,1);1]);
            Poly=[Poly TT1(:).*TT2(:).*TT3(:)];
        end
    end
end

InvCheb=(Poly'*Poly)\Poly';

%% Mapping for Histogram

Gamma = zeros(grid.nb+grid.na+grid.nse,grid.nb+grid.na+grid.nse-3);

for i=1:grid.nb-1
    Gamma(1:grid.nb,i) = - Xss(1:grid.nb);
    Gamma(i,i)         = 1-Xss(i);
    Gamma(i,i)         = Gamma(i,i) - sum(Gamma(1:grid.nb,i));
end

jj = grid.nb;

for i=1:grid.na-1
    Gamma(jj+(1:grid.na),jj+i-1) = - Xss(jj+(1:grid.na));
    Gamma(jj+i,jj-1+i)           = 1-Xss(jj+i);
    Gamma(jj+i,jj-1+i)           = Gamma(jj+i,jj-1+i) - sum(Gamma(jj+(1:grid.na),jj-1+i));
end

jj = grid.nb + grid.na;

for i=1:grid.nse-1
    Gamma(jj+(1:grid.nse),jj+i-2) = - Xss(jj+(1:grid.nse));
    Gamma(jj+i,jj-2+i)            = 1-Xss(jj+i);
    Gamma(jj+i,jj-2+i)            = Gamma(jj+i,jj-2+i) - sum(Gamma(jj+(1:grid.nse),jj-2+i));
end

end
