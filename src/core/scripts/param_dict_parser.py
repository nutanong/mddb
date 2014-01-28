from optparse import OptionParser

def parse(option_dict):
  option_parser = OptionParser()

  for opt_name, default_val in option_dict.iteritems():
    try:
      default_val,comment = default_val
    except:
      comment = opt_name
      pass

    option_parser.add_option(
      "--"+opt_name,
      type=type(default_val), dest=opt_name,
      default=default_val,
      help=comment, metavar="#"+opt_name.upper())

  options,args = option_parser.parse_args()
  return dict(vars(options).items())

