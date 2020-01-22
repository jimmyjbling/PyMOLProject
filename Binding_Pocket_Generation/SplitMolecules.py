import pymol

def split_object(target_obj=None, source_obj=None, max_iter=500, quiet=1, _self=cmd):
    """
DESCRIPTION

    Splits a multi-molecular object into one multi-state object

ARGUMENTS

    target_obj

        (string) name of target object
        
    source_obj

        (string) name of source object
    
    max_iter

        (int) maximum number of object to process; set to 0 to unlimit
    
    """
    if source_object==None:
        print("Error: Please provide a source object.")
        return

    # ensure the user gave us one object; save for prefix

    obj_list = _self.get_object_list(target_obj)

    if len(obj_list)>1:
        print("Error: Please provide only one object at a time.")
        return

    if target_object==None:
        target_object = _self.get_unused_name(source_obj, alwaysnumber=0)

    # grab unused selection name
        
    s = cmd.get_unused_name("_store")

    # make original selection which we'll pare down

    cmd.select(s, source_obj)

    count = 0

    while cmd.count_atoms(s) and count<max_iter:
        count+=1

        # create the object from the first molecular
        # object inside pfx
        cmd.create(pfx, "bm. first " + s, 1, count)

        # remove the first molecular object from
        # the source selection and re-iterate
        cmd.select(s, "%s and not bm. first %s" % (s,s))

    if not quiet:
        print ("Created new object %s." % target_obj)

cmd.extend("split_object", split_object)