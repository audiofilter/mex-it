d = 0.1;
f = single(d);
i8  = int8(-1);
i16 = int16(-2);
i32 = int32(-4);
i64 = int64(-8);
u8 = uint8(1);
u16 = uint16(2);
u32 = uint32(4);
u64 = uint64(8);
b = true;

sum= types_example(d,f,i8,i16,i32,i64,u8,u16,u32,u64,b);

% sum will not be exactly 0.2 due to double/float rounding!
diff = sum - (d + f);

if (diff == 0),
  disp('success');
else
  error('sum is wrong in test_types');
end

% Below are checks that we have the right input types

% Check if we can catch the error
try 
  sum= types_example(d,f,single(i8),i16,i32,i64,u8,u16,u32,u64,b);
catch ME
  switch ME.identifier
    case 'mex_function:validate_and_populate_arg'
      disp ([' caught expected type validation error in input argument with error message: ']);
      disp ([' ',ME.message]);
    otherwise
      error('Unknown error');
  end
end

try 
  sum= types_example(d,f,i8,i16,single(i32),i64,u8,u16,u32,u64,b);
catch ME
  switch ME.identifier
    case 'mex_function:validate_and_populate_arg'
      disp ([' caught expected type validation error in input argument with error message: ']);
      disp ([' ',ME.message]);
    otherwise
      error('Unknown error');
  end
end

try 
  sum= types_example(d,f,i8,i16,i32,i64,u8,single(u16),u32,u64,b);
catch ME
  switch ME.identifier
    case 'mex_function:validate_and_populate_arg'
      disp ([' caught expected type validation error in input argument with error message: ']);
      disp ([' ',ME.message]);
    otherwise
      error('Unknown error');
  end
end

try 
  sum= types_example(d,f,i8,i16,i32,i64,u8,u16,u32,u64,0.1);
catch ME
  switch ME.identifier
    case 'mex_function:validate_and_populate_arg'
      disp ([' caught expected type validation error in input argument with error message: ']);
      disp ([' ',ME.message]);
    otherwise
      error('Unknown error');
  end
end

try 
  sum= types_example(single(d),f,i8,i16,i32,i64,u8,u16,u32,u64,b);
catch ME
  switch ME.identifier
    case 'mex_function:validate_and_populate_arg'
      disp ([' caught expected type validation error in input argument with error message: ']);
      disp ([' ',ME.message]);
    otherwise
      error('Unknown error');
  end
end

