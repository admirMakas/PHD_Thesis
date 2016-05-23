function [bestvar,bestobj,history,eval_count] = ...
    ga(ObjFun,nvars,A,b,Aeq,beq,LB,UB,NONLCON,GAOptions)

% Genetic Algorithm search
% Inputs:
%	ObjFun - character string of function name to be searched
%	nvars - scalar number of variables
%   LB,UB - 1 x k vectors of upper and lower bounds
%   Inputs 3-6 & 9 are not used
% Optional Inputs:
%	GAOptions.PopulationSize - scalar population size (default nvars*10)
%	GAOptions.Generations - scalar number of geneartions (default 50)

%
% Outputs:
%	bestvar - 1 x nvars vector of optimum variables
%	bestobj - scalar optimum of Objfun
%	history - bestobj after each generation
%	eval_count - number of evaluations of Objfun
%
% Copyright 2007 A Sobester and A I J Forrester
%
% This program is free software: you can redistribute it and/or modify  it
% under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or any
% later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
% 
% You should have received a copy of the GNU General Public License and GNU
% Lesser General Public License along with this program. If not, see
% <http://www.gnu.org/licenses/>.

disp('Starting the "Engineering Design via Surrogate Modelling" GA.')

global RANGE_lower RANGE_upper
global VARS
global encoding
global POP

if exist('GAOptions','var')==0
	GAOptions=[];
end

% === RUNTIME PARAMETERS ====

% = Population size =
if isfield(GAOptions,'PopulationSize')==0
	% Default value
    GAOptions.PopulationSize = 50;
end
Popsize = GAOptions.PopulationSize;

% = Number of generations =
if isfield(GAOptions,'Generations')==0
    % Default value
	GAOptions.Generations = 100;
end
Gens = GAOptions.Generations;

% = Tournament size =
if isfield(GAOptions,'TSize')==0
    % Default value
	GAOptions.TSize = 5;
end
Tournamentsize = GAOptions.TSize;

% = Minimization or maximization? =
if isfield(GAOptions,'Direction')==0
    % Default value: minimization
	GAOptions.Direction = -1;
end
Direction = GAOptions.Direction;

% = Encoding - number of bits per variable = 
encoding = 20;

% = Probabilities =
Prep = 0.1;
Pcr = 0.5;
Pmut = 0.4;
Pmutdash = 0.1;

RANGE_lower = LB;
RANGE_upper = UB;

VARS = length(RANGE_lower);

history = [];
ev_hist =[];

%Number of mutated bits
nbits = round(Pmutdash*encoding/Pmut);
eval_count = 0;
failed_count = 0;


% Initial population
if isfield(GAOptions,'POP')==0
	%disp('Generating and evaluating the initial population...')
	j = 0;
	while j < Popsize
		candidate = rand_bin(VARS*encoding);
        cand_fit = ffitness(candidate,ObjFun);
		eval_count = eval_count + 1;
		if ~isnan(cand_fit)
			j = j + 1;
			POPfit(j)=cand_fit;
			POP{j} = candidate;
		else
			failed_count = failed_count + 1;
		end
	end

	init_pop_failed = failed_count;

	ev_hist = [ev_hist,eval_count];
else
	disp('User supplied intial population.')
	for j=1:Popsize
		POP{j}=encode(GAOptions.POP(j,:),RANGE_lower,RANGE_upper,ones(nvars,1).*encoding);
	end
	if isfield(GAOptions,'POPfit')==0
		j = 0;
		InitialPOP=POP;
		while j < Popsize
			if j+failed_count+1<=Popsize	
				candidate = InitialPOP{j+failed_count+1};
				cand_fit = ffitness(candidate,ObjFun);
			else
				candidate = rand_bin(VARS*encoding);
				cand_fit = ffitness(candidate,ObjFun);
			end
			eval_count=eval_count+1;
			if ~isnan(cand_fit)
				j = j + 1
				POPfit(j)=cand_fit;
				POP{j} = candidate;
			else
				failed_count=failed_count+1
			end
		end
	else
		eval_count=0;
	end
	ev_hist=[ev_hist, eval_count];
end

if Direction==-1
    history = [history,min(POPfit)];
else
    history = [history,max(POPfit)];
end



% Main loop iterating over number of generations
for gen = 1:Gens

    if isfield(GAOptions,'target')==1
        if (Direction < 0 && min(POPfit) < GAOptions.target) ||...
           (Direction > 0 && max(POPfit) > GAOptions.target)
           break;
        end
    end
        
    %disp(['Generating and evaluating population ',num2str(gen)])
	
    NEWPOP = {};
	NEWPOPfit = [];
	
    % Make sure best survives
	if Direction < 0
		[mini,ind] = min(POPfit);
	else
		[maxi,ind] = max(POPfit);
    end

    NEWPOP = add_cell(NEWPOP,POP{ind});
	NEWPOPfit = [NEWPOPfit,POPfit(ind)];
    bestsaved = NEWPOPfit(1);
    
    while length(NEWPOP)<Popsize
        % Decide which operator to use
        operator = multiflip(Prep,Pcr);
        switch operator
            case 1
                % Reproduction
                ind = tournament(POPfit,Tournamentsize,Direction);
                NEWPOP = add_cell(NEWPOP,POP{ind});
                NEWPOPfit = [NEWPOPfit,POPfit(ind)];
            case 2
                % Crossover
                ind1 = tournament(POPfit,Tournamentsize,Direction);
                ind2 = tournament(POPfit,Tournamentsize,Direction);
                [off1, off2] = crossover(POP{ind1},POP{ind2},1);
                
                obj1 = ffitness(off1,ObjFun);
            	obj2 = ffitness(off2,ObjFun);

                % Test for evaluation failure
                if isnan(obj1)
                    NEWPOP = add_cell(NEWPOP,POP{ind1});
                    NEWPOPfit = [NEWPOPfit,POPfit(ind1)];
                    failed_count = failed_count + 1;
                else
                        NEWPOP = add_cell(NEWPOP,off1);
                        NEWPOPfit = [NEWPOPfit,obj1];
                end
                        
                if isnan(obj2)
                    NEWPOP = add_cell(NEWPOP,POP{ind2});
                    NEWPOPfit = [NEWPOPfit,POPfit(ind2)];
                    failed_count = failed_count + 1;
                else
                        NEWPOP = add_cell(NEWPOP,off2);
                        NEWPOPfit = [NEWPOPfit,obj2];
                end
                
                eval_count = eval_count + 2;
                
         	case 3
                % Mutation
                ind = tournament(POPfit,Tournamentsize,Direction);
                off = mutate(POP{ind},nbits);
                obj = ffitness(off,ObjFun);
                
                % Test for evaluation failure
	    		if isnan(obj)
                    NEWPOP = add_cell(NEWPOP,POP{ind});
                    NEWPOPfit = [NEWPOPfit,POPfit(ind)];
                    failed_count = failed_count + 1;
                else
                        NEWPOP = add_cell(NEWPOP,off);
                        NEWPOPfit = [NEWPOPfit,obj];
                end
                
                eval_count = eval_count + 1;
        end %switch
    
    end %generating new population
      
   
    % Should there be extra individuals (one or two)...
	if length(NEWPOP)>Popsize
		NEWPOP = NEWPOP(1:Popsize);
		NEWPOPfit = NEWPOPfit(1:Popsize);
	end      
	POP = NEWPOP;
	POPfit = NEWPOPfit;
	ev_hist = [ev_hist,eval_count];
    if Direction < 0
        history = [history,min(POPfit)];
    else
        history = [history,max(POPfit)];
    end
end

for i=1:length(RANGE_lower)
   GL(i) = encoding;
end

ind     = tournament(POPfit,Popsize,Direction);
bestvar = decode(POP{ind},RANGE_lower,RANGE_upper,GL);
bestobj = POPfit(ind);


disp('GA run complete.')

return



function rstr=rand_bin(length)
% RAND_BIN returns a random binary string of specified length
% RAND_BIN(L) returns a random binary string of length L.
rstr='';
rn=rand(1,length);
for poz=1:length
    if rn(poz)<0.5
        rstr=[rstr,'0'];
    else
        rstr=[rstr,'1'];
    end
end 
return



function f = ffitness(bitstring,ObjFun)
global RANGE_lower RANGE_upper
global encoding
% This looks unnecesarily complicated, but it was needed
% for a multi-species implementation

for i=1:length(RANGE_lower)
   GL(i) = encoding;
end

x = decode(bitstring,RANGE_lower,RANGE_upper,GL);
f = feval(ObjFun,x);

return



function ind = tournament(POPfit,Tournamentsize,Direction)
% Direction: 1 = maximization
%            -1 = minimization
pool = [];
for i=1:Tournamentsize
   pool = [pool,int_rnd_no(length(POPfit))];
end
if Direction==1
	[a,b] = max(POPfit(pool));
else
	[a,b] = min(POPfit(pool));
end
ind = pool(b);
return




function [o1, o2] = crossover(p1,p2,type)
% CROSSOVER(P1,P2,type) crossover operator
% p1, p2: parents
% 
switch type
case 1
   % Simple one-point crossover
   cross_point = int_rnd_no(length(p1)-1);
   o1 = [p1(1:cross_point),p2(cross_point+1:length(p2))];
   o2 = [p2(1:cross_point),p1(cross_point+1:length(p2))];
case 2
   % Two-point crossover
   cross_point1 = int_rnd_no(length(p1)-2);
   cross_point2 = int_rnd_no(length(p1)-1);
   if cross_point1==cross_point2
      o1 = [p1(1:cross_point1),p2(cross_point1+1:length(p2))];
      o2 = [p2(1:cross_point1),p1(cross_point1+1:length(p2))];
   else
      if cross_point1 > cross_point2
         [cross_point2, cross_point1] = swap(cross_point2,cross_point1);
      end
      o1 = [p1(1:cross_point1),p2(cross_point1+1:cross_point2),p1(cross_point2+1:length(p2))];
      o2 = [p2(1:cross_point1),p1(cross_point1+1:cross_point2),p2(cross_point2+1:length(p2))];
   end 
end
return




function ca2 = add_cell(ca1, c)
% Adds cell c to cell array ca1

if length(ca1)==1 && isempty(ca1{1})
   ca2 = ca1;
   ca2{1} = c;
else
   ca2 = cell(1,length(ca1)+1);
   ca2(1:length(ca1)) = ca1;
   ca2{length(ca1)+1} = c;
end
return



function indm=mutate(chromosome,bits)
% MUTATE mutates a desired number of bits in a chromosome
% MUTATE(C,B) mutates a number B of bits in the string C
% Note: the mutation loci are randomly generated and it may 
% happen that not all loci are distinct.

l=length(chromosome);
positions=floor(rand(bits,1)*l)+1;
for i=1:bits
    if chromosome(positions(i))=='1'
        chromosome(positions(i))='0';
    else chromosome(positions(i))='1';
    end
end
indm=chromosome;
return



function x=decode(chromosome,LLIMIT,ULIMIT,GL)
% DECODE decodes a binary chromosome
% DECODE(C,L,U,G) is a vector with length(L) elements, containing
% the genes decoded from the chromosome binary string C, given the
% lower limits (array L) and upper limits (array U) of the respective
% variables. The array GL contains the lengths of the genes of the 
% variables. All three input arrays have the same length, equal to 
% the number of variables.

for i=1:length(LLIMIT)
    x(i)=bin2dec(  chromosome( sum(GL(1:i-1))+1 : sum(GL(1:i)) ) );
    x(i)=LLIMIT(i)+x(i)*(ULIMIT(i)-LLIMIT(i))/(2^GL(i)-1);
end
return



function operator = multiflip(Prep,Pcr)
% Stochastic selection of operator to be used
a = rand;
if a < Prep
   operator = 1;
elseif a >= Prep && a < Prep+Pcr
   operator = 2;
else
   operator = 3;
end
return



function wrn = int_rnd_no(n)
% Generates an integer random number between 1 and n
wrn = floor(rand*n)+1;
return


function chromosome = encode(x,LLIMIT,ULIMIT,GL)
%Encodes a vector into a binary chromosome

for i=1:length(x)
  if x(i)>ULIMIT(i)-0.01
	x(i) = ULIMIT(i)-0.01;
  end
end

chromosome = [];

for i=1:length(LLIMIT)
   unit = (ULIMIT(i)-LLIMIT(i))/(2^GL(i));
   xx = round((x(i)-LLIMIT(i))/unit);
   chromosome = [chromosome,dec2bin(xx,GL(i))];   
end

return
