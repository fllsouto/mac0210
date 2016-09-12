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


function [x_r, x_i, it] = newton(f, f_line, x_0)
  delta = 1.e-8;
  x_k = Inf + i;
  x_k_1 = x_0;
  it = 0;
  mx_it = 15;
  abs_x = abs(x_k_1 - x_k);
  m_abs_x = abs_x;

  while (abs_x > delta) && (it <= mx_it)
    % pause ();
    it++;
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
    % TODO: think about the correctness of the attribution performed in the following else block.
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
    % printf("Iter %d!!!\n", it);

  endwhile

  if(it > mx_it)
    x_r = Inf;
    x_i = Inf;
  endif

endfunction

function[] = write_output(m_result)

  color = 1;
  actual_root = m_result(1, 3);
  fact = 64;
  for i = 1:rows(m_result)
    if(actual_root != m_result(i, 3))
      actual_root = m_result(i, 3);
      color++;
    endif

    m_result(i, 3) = mod(fact, color);
    % printf("%f %f %d ", real(m_result(i, 1)), real(m_result(i, 2)), color);
    % disp(m_result(i, 3));
    % printf("\n")

    % printf("%d %d %d\n", real(m_result(i, 1)), real(m_result(i, 2)), color);

  endfor

  m_result = sortrows(m_result, [1, 2]);

  for i = 1:rows(m_result)
    printf("%f %f %d\n", real(m_result(i, 1)), real(m_result(i, 2)), m_result(i, 3));
  endfor

endfunction

function [] = newton_basins(f_x, n)
  step = 0.02;
  C_line = -n:step:(n - step);
  f_x_line = polyder(f_x);
  delta = 1.e-8;
  m_result = [];

  printf("Total of points : %d\n", columns(C_line)*columns(C_line));

  for x = C_line
    for y = C_line
      % printf("i = %d j = %d\n", x, y);
      [x_r, x_i, it] = newton(f_x, f_x_line, (x + y*i));

      if(x_r == Inf || x_i == Inf)
        m_result = [m_result ; [x, y, Inf]];
      elseif(abs(x_r) < delta)
        m_result = [m_result ; [x, y, (x_i*i)]];
      elseif(abs(x_i) < delta)
        m_result = [m_result ; [x, y, x_r]];
      else
        m_result = [m_result ; [x, y, (x_r + x_i*i)]];
      endif

    endfor
  endfor

  % printf("Matrix : \n")
  % disp(sortrows(m_result, [3]))
  write_output(sortrows(m_result, [3]));
  printf("\nTotal of lines: %d", rows(m_result));
endfunction


f_x = [1, 0, 0, 0, -1];
n = 3;
newton_basins(f_x, n);
printf("\n");
