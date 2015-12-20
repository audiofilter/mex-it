a = 1:5;
b = -0.5*a;

% result should be 5 0s
sum= vector_example(a,b);
sum= vector_example(sum,b);
if (any(sum)),
  error 'problem with eigen sum'
else
  disp('success');
end

% Make 2nd argument be the wrong data type
sb = single(b);

% Check if we can catch the error
try 
  sum = vector_example(a,sb);
catch ME
  switch ME.identifier
    case 'mex_function:validate_and_populate_arg'
      disp ('caught expected type validation error in input argument for vector_example(a,b)');
    otherwise
      error('Unknown error');
  end
end
