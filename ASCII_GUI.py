
def complex_selector(d):
	c = 'a'
	res = []
	while True:
		Announcement("Choose measurements to analyze")
		print_table_from_dict(d, leftcol="Option",rightcol="Measurement")
		print("If you're done selecting: press ENTER")
		print(r'To select everything, type: "*"')
		print(r'To select a range of measurements, type "#-#"')
		c = input("Input Mesaurement: ")
		if c == "":
			break
		else:
			if '-' in c:
				[l,r] = c.split('-')
				tr = [i for i in range(int(l),int(r)+1)]
				[d.pop(i) for i in tr]
				res += tr
				continue
			elif "*" in c:
				return list(d.keys())+res
			else:
				d.pop(int(c))
				res.append(c)
	return res
			
		


def print_table_from_dict(d, **kwargs):
    """Takes a dict and turns it into a nice table"""
    """
    ###############################################
    #        name           #         value       #
    ###############################################
    |       item 1          |      98217340192    |
    |=======================+=====================|
    ###############################################
    """
    leftcol = kwargs.pop("leftcol", "Measurement")
    rightcol = kwargs.pop("rightcol", "Value")
    # Example Table

    # Keys of the dictionary are the "names"
    keys = d.keys()
    len_keys = [len(list(str(key))) for key in keys] + [len(list(str(leftcol)))]
    maxlen_keys = max(len_keys) + 2

    # Values are the things to be printed.
    vals = [d[k] for k in keys]
    len_vals = [len(list(str(v))) for v in vals] + [len(list(str(rightcol)))]
    maxlen_vals = max(len_vals) + 2
    top_border = maxlen_vals + maxlen_keys + 3

    # Print header
    print("#" * top_border)

    print(
        ("#{!s:^" + str(maxlen_keys) + "s}#{!s:^" + str(maxlen_vals) + "s}#").format(
            leftcol, rightcol
        )
    )
    print("#" * top_border)
    for _, value in enumerate(keys):
        print(
            (
                "|{!s:^" + str(maxlen_keys) + "s}+{!s:^" + str(maxlen_vals) + "s}|"
            ).format(value, d[value])
        )
        print("|" + "=" * maxlen_keys + "|" + "=" * maxlen_vals + "|")
    print("#" * top_border)

def select_simple(dic):
	choices = ''
	while True:
		print_table_from_dict(dic)
		try:
			choices = int(input("Enter option number: "))
			print("Your choice was", dic[choices], "returning...")
			break
		except KeyboardInterrupt as e:
			print("keyboard inturrupt recived. Breaking.")
			raise KeyboardInterrupt from e
		except ValueError:
			print("Invalid input. Try again.", '\n'*2)
			continue
		except KeyError:
			print("Incorrect key. Make sure option matches that which"+
			      "is in the table.")
			continue
	return dic[choices]


# Facilitates the ASCIIGUI
def dict_selector(dic):
    """An ASCII gui that forces a user to choose an item within a dict"""
    """
    From an asciigui backend of a project I authored a long time ago.
    https://github.com/Tristan-Anderson/Slifer_Lab_NMR_Toolsuite;

    expects a {key0:tupple0, key1:tupple2, (...) }
        tuple = [String, Function]
    such that the String describes the properties Function.

    the "key" is enumerated, and printed in the first column of the output table.
        the key is what's returned from the selector.

    the "String" (first element of tuple) is used as a human-describer for 
        the second element of the tuple

    the "Function" is the programmatic object that the user "chooses".

    """

    def table_formatter(keys, strs):
        """Formats a table that the user selects from"""
        len_keys = [len(key) for key in keys]
        maxlen_keys = max(len_keys) + 2

        keydescribers = strs
        len_keydescribers = [len(k) for k in keydescribers]
        maxlen_key_descrbers = max(len_keydescribers) + 2
        maxlen_key_descrbers = 10 if maxlen_key_descrbers < 10 else maxlen_key_descrbers
        xbar = maxlen_keys + maxlen_key_descrbers + 4 + 8 + 1

        print("#" * xbar)
        print(
            str(
                "{0:1}{1:^9}{0:1}{2:^"
                + str(maxlen_keys)
                + "}{0:1}{3:^"
                + str(maxlen_key_descrbers)
                + "}{0:1}"
            ).format("#", "option#", "name", "describers")
        )
        print("#" * xbar)
        for index, value in enumerate(keys):
            print(
                str(
                    "{0:1}{1:^9}{0:1}{2:^"
                    + str(maxlen_keys)
                    + "}{0:1}{3:^"
                    + str(maxlen_key_descrbers)
                    + "}{0:1}"
                ).format("#", index, value, keydescribers[index])
            )
        print("#" * xbar)
        return True

    # Take the keys of the dictionary the user has given to us
    keys = dic.keys()
    # Wrap the dictionary keywise into int-keys.
    #  {0:{key0:tupple0}, 1:{key1:tupple2}, (...) }
    options = list(range(len(keys)))
    reconsile = dict(zip(options, keys))

    # Fetch the describers
    keydescribers = []
    for k in keys:
        value = dic[k]
        if isinstance(value,list):
            keydescribers.append(dic[k][0])
        else:
            keydescribers.append(dic[k])
    # Print the table of options
    table_formatter(keys, keydescribers)

    # Force the user to choose an enumerated option
    #   Exploting the try-except KeyError
    choices = ""
    while True:
        try:
            # User inputs
            choices = int(input("Enter option number: "))
            print("Your choice was", reconsile[choices], "returning...")
            break
        except KeyboardInterrupt as e:
            print("keyboard inturrupt recived. Breaking.")
            raise KeyboardInterrupt from e
        except ValueError:  # If the user enters something un-intable
            print("Invalid input. Try again.", "\n" * 2)
            continue
        except KeyError:  # If the user enters a number that isn't on the list.
            print(
                "Incorrect Key. Make sure option number\
                 matches that which is in the table of options."
            )
            continue
    return reconsile[choices]

def Announcement(s: str):
    """
    Puts # around s to make it appear as
        a header to an announcement
    """
    lstr = len(s) + 4
    border = "#" * lstr + ""
    a = "# " + s + " #"
    print(border, a, border, sep="\n")
