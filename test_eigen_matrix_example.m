a = magic(3);
b = 2*a;

[out1] = eigen_add(a,b);
diff = 3*a - out1;
if (any(diff)),
  error 'mismatch'
else
  disp('success')
end

% Use wrong data type for 2nd argument below
sb = single(b);

% Check if we can catch the error
try 
  sum = eigen_add(a,sb);
catch ME
  switch ME.identifier
    case 'mex_function:validate_and_populate_arg'
      disp ([' caught expected type validation error in input argument for vector_example(a,b) with error message: ',ME.message]);
    otherwise
      error('Unknown error');
  end
end
