#! /usr/bin/octave -qf
# Universidade de Sao Paulo - USP
# Instituto de Matematica e Estatistica
#
# MAC0210 - Laboratório de Métodos Numéricos
# Exercicio Programa I - Método de Newton
# Prof : Ernesto G. Birgin
#
# Fellipe Souto Sampaio - 7990422
# Renan Fichberg - 7991131


function [x_r, x_i] = newton(f, f_prime, x_0)
  delta = 1.e-8;
  x_k = Inf + i;
  x_k_1 = x_0;
  it = 0;
  mx_it = 15;
  abs_x = abs(x_k_1 - x_k);
  m_abs_x = abs_x;

  while (abs_x > delta) && (it <= mx_it)
    it++;
    x_k = x_k_1;
    f_x_k =  polyval(f, x_k);

    f_prime_x_k =  polyval(f_prime, x_k);

    if f_prime_x_k != 0
      x_k_1 = x_k - (f_x_k/f_prime_x_k);
      abs_x = abs(x_k_1 - x_k);

      if(abs_x < m_abs_x)
        x_r = real(x_k_1);
        x_i = imag(x_k_1);
        m_abs_x = abs_x;
      endif

      if(it > mx_it)
        x_r = Inf;
        x_i = Inf;
      endif
    else
      x_r = Inf;
      x_i = Inf;
      break;
    endif
  endwhile
endfunction

function[] = write_output(m_result)
  fd = fopen("output.txt", "w");

  m_result = sortrows(m_result, [1, 2]);
  for i = 1:rows(m_result)
    fprintf(fd, "%.8f %.8f %.8f\n", real(m_result(i, 1)), real(m_result(i, 2)), m_result(i, 3));
  endfor
  fclose(fd);
endfunction

function [] = newton_basins(f_x, n)
  step = 0.05;
  C_line = -n:step:n;
  f_x_prime = polyder(f_x);
  delta = 1.e-8;
  m_result = [];
  association = [0, [Inf, Inf]];
  integer = 1;

  printf("Points : %d\n", columns(C_line) * columns(C_line));
  for x = C_line
    for y = C_line
      [x_r, x_i] = newton(f_x, f_x_prime, (x + y*i));

      if(abs(x_r) == Inf || abs(x_i) == Inf)
        m_result = [m_result; [x, y, 0]];
      elseif(rows(association) == 1)
        association = [association; [integer, x_r, x_i]];
        m_result = [m_result; [x, y, integer++]];
      else
        index = 1;
        while(index <= rows(association))
          arr = association(index, 1:3);
          if(abs(x_r - arr(2)) <= delta && abs(x_i - arr(3)) <= delta)
            m_result = [m_result; [x, y, arr(1)]];
            break;
          elseif(index++ == rows(association))
            association = [association; [integer, x_r, x_i]];
            m_result = [m_result; [x, y, integer++]];
            break;
          endif
        endwhile
      endif
    endfor
  endfor
  write_output(m_result);
  printf("Lines: %d\n", rows(m_result));
endfunction

f_x = [23, 0, 2, 42, -1, -12];
n = 4;
newton_basins(f_x, n);
printf("\n");