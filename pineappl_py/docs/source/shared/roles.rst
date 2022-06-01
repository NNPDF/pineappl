:orphan:

.. define custom roles and other site-wide stuffs


.. ----
   code
   ----

.. role:: py(code)
   :language: python

.. role:: bash(code)
   :language: bash


.. ----------
   references
   ----------

..  .. role-- eq

.. ----------
   custom css
   ----------

.. raw:: html

   <style type="text/css">
     span.underlined {
       text-decoration: underline;
     }
     .red {
       color:red;
     }
     .green {
       color:green;
     }
     .blue {
       color:blue;
     }
     .purple {
       color:purple;
     }
   </style>

.. role:: underlined
   :class: underlined

.. role:: red
   :class: red

.. role:: green
   :class: green

.. role:: blue
   :class: blue

.. role:: purple
   :class: purple

.. :underlined:`test`


.. raw:: html

   <style type="text/css">
     span.eqref:before {
       content: "Eq. (";
     }
     span.eqref:after {
       content: ")";
     }
   </style>

.. role:: eqref
   :class: eqref
