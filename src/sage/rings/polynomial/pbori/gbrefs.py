import gzip
from io import StringIO
import uu
import re
from types import ModuleType
from .PyPolyBoRi import Polynomial
AUTO = "auto"
SINGLE = "single"


# def ref_file_name(f):
#     name=re.sub("data/","",f)
#     name=sub(r"\.py","",name)
#     l=name.split("/")[:-1]


def reencode_blocks(block_str):
    return str(block_str).replace(",", "_")


def parse_blocks(block_str, data):
    if block_str == AUTO:
        return data.block_start_hints
    if block_str == SINGLE:
        return []
    return [int(i) for i in block_str.split(",")]


def load_ref_raw(s):
    s = re.sub("data/", "", s)
    s = re.sub(r"data\.", "", s)
    s = re.sub(r"\.py", "", s)
    s = re.sub(r"\.", "/", s)

    ref_file = "ref/" + s + ".ref"
    res_f = open(ref_file)
    res = res_f.read()
    res_f.close()
    return res


def load_ref(s, ordering="lp", blocks=SINGLE):
    return load_ref_gz_uu(s, ordering, blocks)


def ordering_suffix(o, blocks=None):
    if o == "lp":
        return ""
    else:
        if re.match("block", o):
            return "." + o + "_" + reencode_blocks(blocks)
        else:
            return "." + o


def number_of_declared_vars(data):
    try:
        return data.number_of_declared_vars
    except AttributeError:
        return data.r.ngens()


def load_ref_gz_uu(s, o, b):
    s = re.sub("data/", "", s)
    s = re.sub(r"data\.", "", s)
    s = re.sub(r"\.py", "", s)
    s = re.sub(r"\.", "/", s)

    ref_file = "ref/" + s + ordering_suffix(o, b) + ".ref.gz.uu"
    res = StringIO()
    uu.decode(ref_file, res)
    res = res.getvalue()
    res = StringIO(res)
    res = gzip.GzipFile(fileobj=res, mode="r").read()
    res = res.replace(" ", "")
    return res


def convert_refs(ref_file_orig):
    content = open(ref_file_orig).read()
    buf_out = StringIO()
    zipped = gzip.GzipFile(filename=ref_file_orig, mode="w", fileobj=buf_out)
    zipped.write(content)
    zipped.close()
    val = buf_out.getvalue()
    out = open(ref_file_orig + ".gz.uu", "w")
    uu.encode(out_file=out, in_file=StringIO(val))
    out.close()


def dyn_generate(content, name):
    module = ModuleType(name)
    import_header = """from .PyPolyBoRi import Variable,Monomial, Polynomial, Ring, OrderCode
from itertools import chain
from .blocks import AlternatingBlock,Block,AdderBlock,if_then,HigherOrderBlock,declare_ring as orig_declare_ring,declare_block_scheme,MacroBlock\n
def declare_ring(blocks, context=None):
  if context is None:
    context=globals()
  return orig_declare_ring(blocks,context)
"""
    exec(import_header + content, module.__dict__)
    if hasattr(module, "ideal"):
        module.ideal = [Polynomial(p) for p in module.ideal]
    return module


def clean_data(data):
    for a in dir(data):
        if a != "r":
            delattr(data, a)


def load_data(file_name, base_dir="./"):
    in_file = file_name
    if not re.match("^data", in_file):
        in_file = "data/" + in_file
    in_file = re.sub(r".py$", "", in_file)
    in_file = re.sub(r"\.", "/", in_file)
    in_file = in_file + ".py"
    in_file = open(base_dir + in_file).read()
    return dyn_generate(in_file, "pb_data")


def load_file(file_name):

    in_file = file_name

    in_file = open(in_file).read()

    return dyn_generate(in_file, "pb_data")
