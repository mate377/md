classdef md
    properties
        val
        var
    end
    properties(Dependent=true)
        std
    end
    methods(Static=true, Access=private, Hidden=true)
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
    methods(Static=true, Access=protected, Hidden=true)
        function out=arr2struct(arr)
            l=length(arr);
            out(l)=struct('val',[]);
            for i=1:l
            out(i).val=arr(i);
            end
        end
    end
    methods(Static=true, Access=public)
        function out = exprInc( expr, vec )
            args=argnames(expr);
            l=length(args);
            assert(l==length(vec),'md:wrongInput:dimension','vec must have same length of args for expr');
            D=sym('D',[1 l]);
            g2=gradient(expr).^2;
            varExpr=symfun(sum(g2.*D.'),[args D]);
            tmp=md.arr2struct([[vec.val] [vec.var]]);
            out=md(eval(expr(vec.val)),eval(varExpr(tmp.val)),'V');
        end
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
     %   function out=struct2md(a)
      %      
       % end
    end
    methods
        function out=get.std(a)
            out=sqrt(a.var);
        end
        function out=md(val,std,mode)
            if (nargin==0)
                out.val=[];
                out.var=[];
                return
            elseif (nargin==1)
                if isa(val,'md')
                    out=val;
                else
                    out=md(val,0,'V');
                end
            elseif (nargin==3)
                if mode=='V'
                    assert((isnumeric(val)&&isnumeric(std)),'md:wrongInput:type','Numeric values for arguments plz')
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
        % overloading
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
        function out=ne(a,b)
            a=md(a); b=md(b);
            c=a-b;
            sigma=1.5;
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
        
        % math functions
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
        
        % other
        function disp(a)
            if isscalar(a)
                fprintf('\t%s : %g','val',a.val);
                fprintf('\n\t%s : %g\n\n','std',a.std);
                %builtin('disp',a);
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
    end
end