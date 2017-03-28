% The class implement a new type of data, md representing a pair of data
% and uncertainty. Perform calculation on the data require propagate
% uncertainty. Keep the data and his uncertainty separate is counter
% logical, error prone and messes up the workspace.
% md has two properties: val and var which are value and variance of the
% data. A dependent property std=(sqrt(var)) is calculated when needed. I
% usually use std, but all calculation involve var, so I prefer to store it
% and retrieve std when needed. The constructor can be called with std or
% var. I've overloaded a lot (not all yet) of operators and functions in
% order to work with md object as standard numerical value. The propagation
% is made using the formula: var(f(x))= sum_i (df/dxi)^2 * var(xi).
% A new object is made using the constructor md. Say that you call 
% md(val,std), if std is a scalar var can be scalar, vector or matrix and
% the resulting object has the dimensions of val, with the same std and var
% for every element. You can have different uncertainties providing for
% example a vector of values and a vector of uncertainties. works also for
% matrix. Obviously the dimensions must agree. You can construct an md
% object starting form a struct with same structure, through the function
% struct2md.
% other ways to call md are through the addition of a third text argument:
% if it is non present the std is handled as std, if you would call 
% md(val,std,'V') the std will be inserted direclty as var. If you call
% md(val,std,'r') std is now an error relative to the value, for example
% with std 0.1 and 'v' option,the error is 0.1*val. Is useful sometimes
% with instruments which give relative errors.
% There is a function which propagate uncertainty on a general
% expression which you should use when you know in advance the physics law
% on which you'll work. Is called exprInc and want a symfun object with the
% law and a vector of md objects from which substitute the values. I've
% overloaded also disp in order to have a nice presentation of md objects.
% by now the code is extendend when new functionality are needed. 

% TODO:
%   - read the MException class doc in order to concatenate exception in stack
%   - try to do a getter for also val and var, which don't require the
%   square brackets to obtain the value
%   - comment on the latest functions created.
%   - work on test script
%   - add regression, chi2 test
%   - function that pretty print a struct of md data. for example, first
%   the scalar parameters and then the vectors. If the vectors are related
%   it prints the vectors side by side. It will be a static method (struct
%   in input) and maybe it ignore the other type of data (non numeric or 
%   md). It could have a a field like "format" that specifies when vectors
%   are related and how. Maybe I can obtain the behavior subclassing struct
%   object.
%   - I need an object containing linear regression parameters so usual
%   operation like find intersections and print lines become short and
%   natural.

classdef md
    properties
        val
        var
    end
    properties(Dependent=true)
        std
    end
    methods(Static=true, Access=private, Hidden=true)
        % I use it in the constructor, it just check the dimensions
        function out=classdim(a)
            if isscalar(a)
                out=1;
            elseif isvector(a) 
                out=2;
            elseif ismatrix(a)
                out=3;
            else
                throw(MException('md:wrongInput:dimension','Only matrix are handled'))
            end
        end
    end
    %methods(Static=true, Access=protected, Hidden=true)
    %end
    methods(Static=true, Access=public)
        % this function takes in input a symfun object and a vector in
        % order to perform the substitution and obtain finally a value. The
        % propagation of uncertainty is done automatically. Called with
        % only one arg it print the symbolic expression for value and
        % variance.
        function out = exprInc( expr, vec )
            args=argnames(expr);
            l=length(args);
            D=sym('D',[1 l]);
            g2=gradient(expr).^2;
            varExpr=symfun(sum(g2.*D.'),[args D]);
            if nargin==1
                fprintf('val = '); disp(expr);
                fprintf('var = '); disp(varExpr);
                return
            end
            assert(l==length(vec),'md:wrongInput:dimension','vec must have same length of args for expr');
            % I need a cell for pass different arguments to the symfun
            % varExpr
            tmpc=num2cell([[vec.val] [vec.var]]);
            out=md(eval(expr(vec.val)),eval(varExpr(tmpc{:})),'V');
        end
        % transform an md object in a struct with the same structure.
        % Useful as interface between functionality.
        function out=md2struct(a)
            assert(isa(a,'md'),'md:wrongInput:type','type md in input plz');
            [m,n]=size(a);
            out(m,n)=struct('val',[],'std',[]);
            for i=1:m
                for j=1:n
                    out(i,j).val=a(i,j).val;
                    out(i,j).std=a(i,j).std;
                end
            end
        end
        % does the opposite of the precedent method
        function out=struct2md(a)
            try
                fm=cell2mat(fields(a));
            catch ME
                if(ME.identifier=='MATLAB:UndefinedFunction')
                    cause=MException('md:wrongInput:type','need a struct in input');
                    ME.addCause(ME,cause);
                end
                rethrow(ME);
            end
            assert(string(fm(1,:))=='val' && string(fm(2,:))=='std','md:wrongInput:type','struct must have val and std parameters');
            [m,n]=size(a);
            out(m,n)=md();
            for i=1:m
                for k=1:n
                    out(i,k)=md(a(i,k).val,a(i,k).std);
                end
            end
        end
        % it call builtin errorbar with val and std, instead of separate
        % manually the args. Use the args variable, obtainde from the
        % conversion cell2struct allows me to obtain a comma separated list.
        % This csl can be passed as argument to errorbar.
        function out=errorbar(x,y,varargin)
           x=md(x); y=md(y);
           args=cell2struct(varargin,'val',1);
           out=errorbar([x.val],[y.val],[y.std]./2,[y.std]./2,[x.std]./2,[x.std]./2,args.val);
        end
    end
    methods
        % getter for std
        function out=get.std(a)
            out=sqrt(a.var);
        end
        % CONSTRUCTOR 
        % the usual ways to call it are described at the beginning. Called
        % with no args returns an md object with val and var set to []. Is
        % useful in preallocating the space. With one arg is the same as a
        % normal invocation with std=0. It is useful when i have to combine
        % md objects with generic numerical values: a number is well
        % represented as a value with no uncertainty on it. All the check
        % on input is done within the handler of the mode 'V' and produce
        % directly the object with val=val and var=var. The standard call
        % md(val,std) result in a call of md(val,std.^2,'V') which is the
        % same.
        function out=md(val,std,mode)
            if (nargin==0)
                % empty object
                out.val=[];
                out.var=[];
                return
            elseif (nargin==1)
                % construct an object with 0 var
                if isa(val,'md')
                    out=val;
                else
                    out=md(val,0,'V');
                end
            elseif (nargin==3)
                % Here the check on the values in input
                if mode=='V'
                    assert((isnumeric(val)&&isnumeric(std)),'md:wrongInput:type','Numeric values for arguments plz')
                    % go look at the code of classdim, is straightforward
                    nv=md.classdim(val); ns=md.classdim(std);
                    if (nv==1 && ns==1)
                        out.val=val;
                        out.var=std;
                    elseif (nv==2 && ns==1)
                        out(size(val))=md;
                        for i=1:length(val)
                            out(i).val=val(i);
                            out(i).var=std;
                        end
                    elseif (nv==2 && ns==2)
                        assert(isempty(find((size(val)==size(std))==0,1)),'md:wrongInput:dimension','val and std have NOT compatible dimensions');
                        out(size(val))=md;
                        for i=1:length(val)
                            out(i).val=val(i);
                            out(i).var=std(i);
                        end
                    elseif (nv==3 && ns==1)
                        [m,n]=size(val);
                        out(m,n)=md;
                        for i=1:m
                            for j=1:n
                                out(i,j).val=val(i,j);
                                out(i,j).var=std;
                            end
                        end
                    elseif (nv==3 && ns==3)
                        assert(~isempty(find(size(val)==size(std),1)),'md:wrongInput:dimension','value and std have NOT compatible dimensions');
                        [m,n]=size(val);
                        out(m,n)=md;
                        for i=1:m
                            for j=1:n
                                out(i,j).val=val(i,j);
                                out(i,j).var=std(i,j);
                            end
                        end
                    else
                        throw(MException('md:wrongInput:dimension','value and std have NOT compatible dimensions'));
                    end
                elseif mode=='r'
                        assert(isempty(find(std<0,1)),'md:wrongInput:value','std must be positive');
                        out=md(val,(val.*std).^2,'V');
                else
                        throw(MException('md:wrongInput:mode','Unknown mode'));
                end
            elseif (nargin>3)
                throw(MException('md:wrongInput:nargin','Too many input arguments'));
            else % nargin ==2
                assert(isempty(find(std<0,1)),'md:wrongInput:value','std must be positive');
                out=md(val,std.^2,'V');
            end
        end
        % OVERLOADING
        % here i overload the operators. Look at the docs for see at what
        % symbol the function correspond. for example '+' calls 'plus'. The
        % pattern is the same for every function: first I pass the objects
        % through the constructor (if the object is already md nothing
        % changes, if is not md is promoted to md) and then i propagate the
        % uncertainty. Not all the operators have been yet overloaded.
        % Maybe someone is not worth it. The association is obvious for all
        % (i mean: the plus is a plus) but for the ne, which I defined as
        % compatibility operator. Go look at the code.
        function out=uplus(a)
            out=md([a.val],[a.var],'V');
        end
        function out=plus(a,b)
            a=md(a); b=md(b);
            out=md([a.val]+[b.val],[a.var]+[b.var],'V');
        end
        function out=minus(a,b)
            a=md(a); b=md(b);
            out=md([a.val]-[b.val],[a.var]+[b.var],'V');
        end
        function out=uminus(a)
            out=md(-[a.val],[a.var],'V');
        end
        function out=times(a,b)
            a=md(a); b=md(b);
            out=md([a.val].*[b.val],[b.val].^2.*[a.var]+[a.val].^2.*[b.var],'V');
        end
        function out=mtimes(a,b)
            a=md(a); b=md(b);
            out=md([a.val]*[b.val],[b.val].^2*[a.var]+[a.val].^2*[b.var],'V');
        end
        function out=rdivide(a,b)
            a=md(a); b=md(b);
            invb2=[b.val].^(-2);
            a2=[a.val].^2;
            out=md([a.val]./[b.val],invb2.*([a.var]+[b.var].*a2.*invb2),'V');
        end
        % ldivide
        function out=mrdivide(a,b)
            a=md(a); b=md(b);
            invb2=[b.val].^(-2);
            a2=[a.val].^2;
            out=md([a.val]/[b.val],invb2.*([a.var]+[b.var].*a2.*invb2),'V');
        end
        % mldivide
        % lt
        % gt
        % le
        % ge
        
        % ne normally represent the 'not equal' operator. Here it represent
        % the compatible operator, that is: it returns one if the data are
        % compatible, otherwise it returns 0. If you want to change the
        % sigma you must come here. Maybe it will be another property that
        % you can set, probably hidden, but i used only here. I dunno.
        function out=ne(a,b)
            a=md(a); b=md(b);
            sigma=1.5; % 1 is too few, 2 is too much, how many princess are
            % there in the castle? No, seriously.
            c=a-b;
            out=abs([c.val])<sigma.*[c.inc];
        end
        % eq
        % and
        % or
        % not
        % colon
        function out=ctranspose(a)
            out=md([a.val]',[a.var]','V');
        end
        function out=transpose(a)
            out=md([a.val].',[a.var].','V');
        end
        % horzcat
        % vertcat
        % subsref
        % subasgn
        % subsindex
        
        % MATH FUNCTIONS
        % well here they are very few. I'm gonna implement them as
        % needed...
        function out=exp(a)
            expa=exp([a.val]);
            out=md(expa,expa.*[a.var],'V');
        end
        function out=sin(a)
            out=md(sin([a.val]),cos([a.val]).*[a.var],'V');
        end
        function out=cos(a)
            out=md(cos([a.val]),-sin([a.val]).*[a.var],'V');
        end
        function out=log(a)
            out=md(log([a.val]),([a.val].^-2).*[a.var],'V');
        end
        % Here i've overload the disp function which is the function
        % automatically called when the interpreter need to display the
        % object. For a scalar is not very different from the default, but
        % for the vectors it print a table: 1column is the progressive
        % number, 2column is the value and 3column the uncertainty. or the
        % matrix it simply print first the matrix of the value and then the
        % matrix of the std.
        function disp(a)
            if isscalar(a)
                fprintf('\t%s : %g','val',a.val);
                fprintf('\n\t%s : %g\n\n','std',a.std);
            elseif isvector(a)
                str=sprintf('  %c   %12s   %12s\n','n','value','std');
                builtin('disp',str);
                l=length(a);
                for i=1:l
                    fprintf('%3d   %+12.5g   %+12.5g\n',i,a(i).val,a(i).std);
                end
                fprintf('\n');
            elseif ismatrix(a)
                disp('value=\n');
                builtin('disp',[a.val]);
                disp('std=\n');
                builtin('disp',[a.std]);
            else
                throw(MException('md:programFlow:unreachablePoint','Something really bad has happened'));
            end
        end
        % Lchi2 perform a chi2 linear regression of data md in input and
        % calculate the parameters of the low in the form A+B*X. The
        % uncertainties on the x are first discarded  and then propagated
        % through the linear law.
        function [A, B, chi] = Lchi2(x,y,n)

            x1 = [x.val];
            dx1 = [x.std];
            y1 = [y.val];
            dy1 = [y.std];

            dyt = dy1;

            for k=1:n
            w = dyt.^(-2);
            tB = (sum(w)*sum(w.*x1.*y1) - sum(w.*y1)*sum(w.*x1))/(sum(w)*sum(w.*(x1.^2)) - (sum(w.*x1))^2);
            dyt = sqrt(dyt.^2 + (tB*dx1).^2);
            end

            tA = (sum(w.*(x1.^2))*sum(w.*y1) - sum(w.*x1)*sum(w.*x1.*y1))/(sum(w)*sum(w.*(x1.^2)) - (sum(w.*x1))^2);
            tdA = sqrt(sum(w.*(x1.^2))/(sum(w)*sum(w.*(x1.^2)) - (sum(w.*x1))^2));
            tdB = sqrt(sum(w)/(sum(w)*sum(w.*(x1.^2)) - (sum(w.*x1))^2));
            chi = sum(((y1 - tA*ones(1, length(y1))-tB*x1).^2)./(dyt.^2));
            A=md(tA,tdA);
            B=md(tB,tdB);
        end
    end
end