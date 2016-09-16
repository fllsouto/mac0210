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

function [] = print_root(sub_interval, sub_lower_edge, sub_higher_edge, root)
  printf("Interval #%d: ]%g, %g] has root %.8f\n", sub_interval, sub_lower_edge, sub_higher_edge, root);
endfunction

function [] = get_root(ev_f_x, ev_f_x_prime, x, choice, tol, step, new_sign, sub_interval, sub_lower_edge, sub_higher_edge)
  got_root = false;

  while(!got_root)
    % Apply the Newton Method
    root = newton(ev_f_x, ev_f_x_prime, x, choice, tol);

    % checks if the method failed or not
    if(root == Inf)

      % Apply the Bisection method to find a better initial value
      x = bisection(new_sign, x, step, choice)

    else
      got_root = true;
      print_root(sub_interval, sub_lower_edge, sub_higher_edge, root);
    endif
  endwhile
endfunction

function [x_mid] = bisection(new_sign, x, step, choice)
  x_mid = 1;

  % **Recall** that we know that the last iteration x have an opposite sign to actual x
  if(new_sign == 1)
    x_with_pos_f_x = x; x_with_neg_f_x = x - step;
  else
    x_with_pos_f_x = x - step; x_with_neg_f_x = x;
  endif

  % Apply the bisection method 3 times in a row.
  for it = 1:1:3

    % get the mean value from the two distinct x values
    x_mid = (x_with_pos_f_x + x_with_neg_f_x) / 2;
    [ev_f_x_mid, ev_f_x_mid_prime, ev_f_x_mid_double_prime] = evaluate_by_choice(choice, x_mid);

    % Update the correct value
    if(ev_f_x_mid >= 0)
      x_with_pos_f_x = x_mid;
    else
      x_with_neg_f_x = x_mid;
    endif
  endfor
endfunction

function [root] = newton(f, f_prime, x_k, choice, tol)
  x_k_1 = x_k - (f / f_prime);
  root = x_k_1;
  [ev_f_x_k_1, ev_f_x_k_1_prime, ev_f_x_k_1_double_prime] = evaluate_by_choice(choice, x_k_1);

  while(abs(ev_f_x_k_1) >= tol)
    x_k = x_k_1;
    ev_f_x_k = ev_f_x_k_1;
    ev_f_x_k_prime = ev_f_x_k_1_prime;
    x_k_1 = x_k - (ev_f_x_k / ev_f_x_k_prime);
    [ev_f_x_k_1, ev_f_x_k_1_prime, ev_f_x_k_1_double_prime] = evaluate_by_choice(choice, x_k_1);
    % compare last evaluations to see if we are getting significant results
    if(abs(ev_f_x_k_1) > abs(ev_f_x_k) / 2)
      root = Inf;
      break;
    else
      root = x_k_1;
    endif
  endwhile
endfunction

function [] = intervals(choice, interval, n_inter, tol)
  % X axis distance between each value of x to be evaluated
  step = 0.05;

  % All values of x to be evaluated and contained in this range.
  % Limits will be verified independently.
  range = (interval(1) + step):step:(interval(2) - step);

  % tolerance
  delta = 1.e-8;

  % sub interval counter where 1 <= sub_interval <= n_inter.
  sub_interval = 1

  % first sub interval
  sub_inter_length = (interval(2) - interval(1)) / n_inter;
  sub_lower_edge = interval(1);
  sub_higher_edge = sub_lower_edge + sub_inter_length;
  % first sub interval is [a, c], second is ]c, d], third ]d, e]
  % and so on till last ]z , b]. **Recall whole interval is [a, b]**

  % get initial sign of the function f_x. 1 == positive, -1 == negative.
  [ev_f_x, ev_f_x_prime, ev_f_x_double_prime] = evaluate_by_choice(choice, interval(1));
  if(ev_f_x > 0)
    sign = 1;
    new_sign = 1;
  elseif(ev_f_x < 0)
    sign = -1;
    new_sign = -1;
  else % Lucky! It is a root
    sign = 1;
    new_sign = 1;
    print_root(sub_interval, sub_lower_edge, sub_higher_edge, interval(1));
  endif

  printf("Points: %d\n", columns(range) + 2); % +2 considering the limits that are not computed in the loop

  for x = range
    % checks if x has entered a new sub interval
    if(x > sub_higher_edge)
      sub_interval++;
      sub_lower_edge = sub_higher_edge;
      sub_higher_edge = sub_lower_edge + sub_inter_length;
    endif

    % Evaluate the correct function on point x
    [ev_f_x, ev_f_x_prime, ev_f_x_double_prime] = evaluate_by_choice(choice, x);

    % Get sign of evaluated f(x) for the x value
    if(ev_f_x > 0)
      new_sign = 1;
    elseif(ev_f_x < 0)
      new_sign = -1;
    else % Lucky! It is a root
      new_sign = 1;
      print_root(sub_interval, sub_lower_edge, sub_higher_edge, x);
      continue;
    endif

    % check if our function has passed through an undiscovered root
    if(new_sign != sign)
      get_root(ev_f_x, ev_f_x_prime, x, choice, tol, step, new_sign, sub_interval, sub_lower_edge, sub_higher_edge);
      sign = new_sign;
    endif
  endfor

  % Checks the last point of the range
  [ev_f_x, ev_f_x_prime, ev_f_x_double_prime] = evaluate_by_choice(choice, interval(2));
  if(ev_f_x > 0)
    new_sign = 1;
  elseif(ev_f_x < 0)
    new_sign = -1;
  else % Lucky! It is a root
    new_sign = 1;
    print_root(sub_interval, sub_lower_edge, sub_higher_edge, interval(2));
  endif

  % Check if limit of the range has a root (that is, the last step between the last x and the one before it).
  if(new_sign != sign)
    get_root(ev_f_x, ev_f_x_prime, interval(2), choice, tol, step, new_sign, sub_interval, sub_lower_edge, sub_higher_edge);
  endif
endfunction

% -------------------------------------------
% Enunciations Math Functions
% -------------------------------------------
% Choice 1: f(x) = 2cosh(x/4) - x

% function 1 from enunciation
function [f_x] = function1(x)
  f_x = 2 * cosh(x / 4);
  f_x = f_x - x;
endfunction

% first derivative of function 1 from enunciation
function [f_x_prime] = function1_der1(x)
  f_x_prime = sinh(x / 4) / 2;
  f_x_prime = f_x_prime - 1;
endfunction

% second derivative of function 1 from enunciation
function [f_x_double_prime] = function1_der2(x)
  f_x_double_prime = cosh(x / 4) / 8;
endfunction

% -------------------------------------------
% Choice 2: f(x) = sin(x) / x if x != 0 or f(x) = 1 if x = 0

% function 2 from enunciation
function [f_x] = function2(x)
  if(x == 0)
    f_x = 1;
  else
    f_x = sin(x) / x;
  endif
endfunction

% first derivative of function 2 from enunciation
function [f_x_prime] = function2_der1(x)
  if(x == 0)
    f_x_prime = 1;
  else
    f_x_prime = ((x * cos(x)) - sin(x)) / (x * x);
  endif
endfunction

% second derivative of function 2 from enunciation
function [f_x_double_prime] = function2_der2(x)
  if(x == 0)
    f_x_double_prime = 1;
  else
    f_x_double_prime = -2 * ((x * cos(x)) - sin(x));
    f_x_double_prime = f_x_double_prime / (x * x * x);
    f_x_double_prime = f_x_double_prime - (sin(x) / x);
  endif
endfunction

% -------------------------------------------

% Polynomial function
function [f_x] = function3(x)
  f_x = [1, 2, -3];
endfunction

% Polynomial function first derivative
function [f_x_prime] = function3_der1(f_x)
  f_x_prime = polyder(f_x);
endfunction

% Polynomial function second derivative
function [f_x_double_prime] = function3_der2(f_x_prime)
  f_x_double_prime = polyder(f_x_prime);
endfunction

% -------------------------------------------
% Mathematical Function getter based on choice.
% Returns evaluated values of the point for the
% function and its first and second derivative

function [ev_f_x, ev_f_x_prime, ev_f_x_double_prime] = evaluate_by_choice(choice, x)
  if(choice == 1)
    ev_f_x = function1(x);
    ev_f_x_prime = function1_der1(x);
    ev_f_x_double_prime = function1_der2(x);
  elseif(choice == 2)
    ev_f_x = function2(x);
    ev_f_x_prime = function2_der1(x);
    ev_f_x_double_prime = function2_der2(x);
  else
    f_x = function3(x);
    f_x_prime = function3_der1(f_x);
    f_x_double_prime = function3_der2(f_x_prime);
    ev_f_x = polyval(f_x, x);
    ev_f_x_prime = polyval(f_x_prime, x);
    ev_f_x_double_prime = polyval(f_x_double_prime, x);
  endif
endfunction

% -------------------------------------------
% Program start

# choice can be either 1, 2 or 3.
# 1 for enunciation function 1
# 2 for enunciation function 2
# 3 for polynomials of degree 2 or greater
choice = 3;

% function interval
lower_edge = -10; % a from enunciation
higher_edge = 10; % b from enunciation
interval = [lower_edge, higher_edge];

% Number of sub intervals
n_inter = 10;

# tolerance to accept calculated value by newton method as root
tol = 1.e-8;

if (choice == 1 || choice == 2 || choice == 3)
  intervals(choice, interval, n_inter, tol);
else
  printf("Unexpected 'choice' value. It must be 1, 2 or 3.\n\n");
  printf("Choice 1: f(x) = 2cosh(x/4) - x\n");
  printf("Choice 2: f(x) = sin(x) / x if x != 0 or f(x) = 1 if x = 0\n");
  printf("Choice 3: f(x) is a polynomial with a degree of 2 or greater\n");
endif
