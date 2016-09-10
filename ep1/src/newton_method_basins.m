#! /bin/octave -qf
#Universidade de Sao Paulo - USP
#Instituto de Matematica e Estatistica
#
#MAC0210 - Laboratório de Métodos Numéricos
#Exercicio Programa I - Método de Newton
#Prof : Ernesto G. Birgin
#
#Fellipe Souto Sampaio - 7990422
#Renan Fichberg




function [x_r, x_i, in] = newton(f, f_line, x_0)
  delta = 1.e-8;
  x_k = Inf + i;
  x_k_1 = x_0;
  in = 0;
  mx_in = 30;
  abs_x = abs(x_k_1 - x_k);
  m_abs_x = abs_x;

  while (abs_x > delta) && (in < mx_in)
    % pause ();
    in++;
    x_k = x_k_1;
    f_x_k =  polyval(f, x_k);
    
    % Divisao por zero em potencial
    f_line_x_k =  polyval(f_line, x_k);


    % printf("f_line_x_k : %f\n")
    % disp( f_line_x_k)

    % printf("f_x_k : %f\n")
    % disp(f_x_k)

    if f_line_x_k != 0
      f_q = f_x_k/f_line_x_k;
    else
      f_q = f_x_k;
    endif

    x_k_1 = x_k - (f_q);
    
    % printf("x_k : %f\n")
    % disp(x_k)

    % printf("x_k_1 : %f\n")
    % disp(x_k_1)

    % printf("abs: %f\n\n\n", abs(x_k_1 - x_k))
    abs_x = abs(x_k_1 - x_k);

    if(abs_x < m_abs_x)
      x_r = real(x_k_1);
      x_i = imag(x_k_1);
      m_abs_x = abs_x;
    endif
    % printf("Inter %d!!!\n", in);

  endwhile
  % printf("After For!!!\n");

endfunction

function [] = newton_basins(f_x, n)
  step = 0.25;
  C_line = -n:step:n;
  f_x_line = polyder(f_x);
  delta = 1.e-8;
  
  for x = C_line
    for y = C_line
      % printf("i = %d j = %d\n", x, y);
      [x_r, x_i, in] = newton(f_x, f_x_line, (x + y*i));
      if(abs(x_r) < delta)
        printf("%d %d %f\n", x, y, x_i);
      else
        printf("%d %d %f\n", x, y, x_r);
      endif
    endfor
  endfor

endfunction


f_x = [1, 0, 0, 0, -1];
n = 2;
newton_basins(f_x, n);