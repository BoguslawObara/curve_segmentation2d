function S = curve_segmentation2d(M) 
%% curve_segmentation2d - linear segmentation of digital curve
%   
%   REFERENCE:
%       I. Debled-Rennesson and J.P. Reveilles,
%       A linear algorithm for segmentation of digital curves,
%       International Journal of Pattern Recognition and Artificial 
%       Intelligence, 9, 635-662,1995
%
%       A. Lichius, A. B. Goryachev, M. D. Fricker, B. Obara, 
%       E. Castro-Longoria, N. D. Read, CDC-42 and RAC-1 regulate opposite 
%       chemotropisms in Neurospora crassa, 
%       Journal of Cell Science, 127, 9, 1953-1965, 2014
%
%   INPUT:
%       M     - contour as M = [y,x].
%
%   OUTPUT:
%       S     - linear segments
%
%   AUTHOR:
%       Boguslaw Obara

%% RUN
x0 = M(1,1); y0=M(1,2);
N = size(M,1);

dM = compute_dM(M);
S = segmentation(x0,y0,dM);

%% SUB-FUNCTIONS
%% plot
    function plot_segments(M,S)
        plot(M(:,1), M(:,2),'*'); xlim([0,max(M(:,1))+1]); ylim([0,max(M(:,2))+1]); hold on
        plot(M(S(:),1), M(S(:),2),'ro','MarkerEdgeColor','r','MarkerSize',10); 
    end
%% dM
    function dM = compute_dM(M)
        for i=1:size(M,1)-1
            dM.x(i) = M(i+1,1) - M(i,1);
            dM.y(i) = M(i+1,2) - M(i,2);
        end
        dM.x = [0, dM.x];
        dM.y = [0, dM.y];
    end
%% segmentation
    function S = segmentation(x0,y0,dM)
        %initialization();
        use_separ_jointing = 0;
        k = 1; xk = 0; yk = 0;  % k=0; -> but Matlab k=1;
        alpha = 0; beta = 0; gamma = 0; delta = 0; phi = 0; psi = 0;
        a = 0; b = 0; mi = 0; Ux = 0; Uy = 0; Lx = 0; Ly = 0;
        h = 1; xh = x0; yh = y0; % h=0; -> but Matlab h=1;
        S = k;
        while h < N
            [k,xk,yk,alpha,beta,gamma,delta,phi,psi,a,b,mi,Ux,Uy,Lx,Ly] = move_frame(k,xk,yk,alpha,beta,gamma,delta,phi,psi,a,b,mi,Ux,Uy,Lx,Ly);
            [k,xk,yk,a,b,mi] = recognize_segment(k,xk,yk,a,b,mi);
            %disp(['Segment: ' sprintf('%f %f %f %f %f', [h k mi a b])])
            S = [S;k]; 
            if use_separ_jointing == 1
                h = k; xh = xk; yh = yk; % joining points
            else
                if k == N
                    h = N;
                else
                    h = k + 1; xh = xk + dM.x(h); yh = yk + dM.y(h);
                end
            end
        end
        %% move_frame
        function [k,xk,yk,alpha,beta,gamma,delta,phi,psi,a,b,mi,Ux,Uy,Lx,Ly] = move_frame(k,xk,yk,alpha,beta,gamma,delta,phi,psi,a,b,mi,Ux,Uy,Lx,Ly);
            k = h; xk = xh; yk = yh;
            dx = dM.x(h+1); dy = dM.y(h+1);
            while (k < N) && (dM.x(k+1) == dx) && (dM.y(k+1) == dy)
                k = k + 1; xk = xk + dM.x(k); yk = yk + dM.y(k);
            end
            [alpha,beta,gamma,delta,phi,psi] = rotation(xk,yk,dx,dy,alpha,beta,gamma,delta,phi,psi);
            if k < N
                [alpha,beta,gamma,delta,phi,psi] = adjust_axis(dx,dy,dM.x(k+1),dM.y(k+1),alpha,beta,gamma,delta,phi,psi);
            end
            if (dx == 0) || (dy == 0)
                Lx = k - h;
            elseif Lx == 0
                a = 1; b = Lx + 1; mi = 0; Ly = 0; Ux = 0; Uy = 0;
            end
        end
        %% rotation
        function [alpha,beta,gamma,delta,phi,psi] = rotation(xh,yh,dx,dy,alpha,beta,gamma,delta,phi,psi)
            if (dx == 1) && (dy >= 0)
                alpha = 1;  beta = 0;  phi = -xh;
                gamma = 0; delta = 0; psi = -yh;
            elseif (dx <= 0) && (dy == 1)
                alpha = 0;  beta = 1;  phi = yh;
                gamma = -1; delta = 0; psi = xh;
            elseif (dx == -1) && (dy <= 0)
                alpha = -1;  beta = 0;  phi = xh;
                gamma = 0; delta = -1; psi = yh;
            elseif (dx >= 0) && (dy == -1)
                alpha = 0;  beta = -1;  phi = yh;
                gamma = 1; delta = 0; psi = -xh;
            end
        end
        %% adjust_axis
        function [alpha,beta,gamma,delta,phi,psi] = adjust_axis(dx,dy,ddx,ddy,alpha,beta,gamma,delta,phi,psi)
            if (dy == 0) && (ddx == dx) && (ddy == -dx)
                [alpha,beta,gamma,delta,phi,psi] = change_vertical_axis(alpha,beta,gamma,delta,phi,psi);
            elseif (dx == 0) && (ddx == dy) && (ddy == dy)
                [alpha,beta,gamma,delta,phi,psi] = change_horizontal_axis(alpha,beta,gamma,delta,phi,psi);
            else
                if ((dx == dy) && (ddx == 0) && (ddy == dy)) || ((dx == -dy) && (ddx == dx) && (ddy == 0))
                    [alpha,beta,gamma,delta,phi,psi] = exchange_axis(alpha,beta,gamma,delta,phi,psi);
                end
            end
        end
        %% change_vertical_axis
        function [alpha,beta,gamma,delta,phi,psi] = change_vertical_axis(alpha,beta,gamma,delta,phi,psi)
            if delta ~= 0
                delta = -delta; psi = -psi;
            elseif beta ~= 0
                beta = -beta; phi = -phi;
            end
        end
        %% change_horizontal_axis
        function [alpha,beta,gamma,delta,phi,psi] = change_horizontal_axis(alpha,beta,gamma,delta,phi,psi)
            if alpha ~= 0
                alpha = -alpha; phi = -phi;
            elseif gamma ~= 0
                gamma = -gamma; psi = -psi;
            end
        end
        %% exchange_axis
        function [alpha,beta,gamma,delta,phi,psi] = exchange_axis(alpha,beta,gamma,delta,phi,psi)
            temp = alpha; alpha = gamma; gamma = temp;
            temp = beta;  beta = delta;  delta = temp;
            temp = phi;   phi = psi;     psi = temp;        
        end
        %% recognize_segment
        function [k,xk,yk,a,b,mi] = recognize_segment(k,xk,yk,a,b,mi)
            while k<N
                ddx = dM.x(k+1); ddy = dM.y(k+1);
                ddx_loc = alpha*ddx + beta*ddy;
                ddy_loc = gamma*ddx + delta*ddy;
                if (ddx_loc <= 0) || (ddy_loc < 0) %% then live - break?
                    break;
                end
                xx = xk + ddx; yy = yk + ddy;
                x_loc = alpha*xx + beta*yy + phi;
                y_loc = gamma*xx + delta*yy + psi;
                r = a*x_loc - b*y_loc;
                if (r < mi-1) || (r > mi+b)
                    break;
                end;
                k = k+1; xk = xx; yk = yy;
                if r == (mi+b)
                    Ux = x_loc - b; Uy = y_loc + 1 - a;
                    a = y_loc - Ly; b = x_loc - Lx;
                    mi = a*Lx - b*Ly - b + 1;
                elseif r == (mi-1)
                    Lx = x_loc - b; Ly = y_loc - 1 - a;
                    a = y_loc - Uy; b = x_loc - Ux;
                    mi = a*Ux - b*Uy;
                end
            end
        end
        %% initialization
        function [] = initialization()
            M = M - 1;
            k = 1;
            while k<N
                l = M(k,:) - M(1,:);
                r = (k-1)*(M(2,:) - M(1,:));
                if r==l
                    k = k + 1;
                else
                    k = k - 1;
                    break;
                end
            end
            dx  =  M(2,1) - M(1,1);
            dy  =  M(2,2) - M(1,2);
            if (dx == 1) && (dy == 0)
                a = 1; b = ((k-1) + 1); mi = 0;%((k-1) + 1);
                U = M(1,:); L = M(k,:);
            elseif (dx == 1) && (dy == 1)
                a = 1; b = 1; mi = 0;
                U = M(1,:); L = M(k,:);                
            end
            %%--------
            xU = U(1); yU = U(2); xL = L(1); yL = L(2);
            while k<N
                dx  =  M(k+1,1) - M(k,1);
                dy  =  M(k+1,2) - M(k,2);
                xM = M(k,1); yM = M(k,2);
                if (dx <= 0) || (dy < 0); break; end
                xx = xM + dx; yy = yM + dy;
                r = a*xx - b*yy;
                if (r < (mi-1)) || (r > (mi+b)); break; end
                k = k + 1; xM = xx; yM = yy;
                if r == (mi+b)
                    xU = xM - b; yU = yM + 1 - a;
                    a = yM - yL; b = xM - xL;
                    mi = a*xM - b*yM - b + 1;
                elseif r == (mi-1)
                    xL = xM - b; yL = yM - 1 - a;
                    a = yM - yU; b = xM - xU;
                    mi = a*xM - b*yM;                
                end
                U = [xU yU]; L = [xL yL];
                %M;
            end
        end
    %% segmentation end
    end
%%
end