def get_3d():

    a = """<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/1.8.0/3Dmol.js"></script>

    <style>
        .mol-container {
        width: 100%;
        height: 100%;
        position: relative;
        }
    </style>

    <div id="container-01" class="mol-container"></div>

    <script>
        $(function() {
        
        let element = $('#container-01');
        let viewer = $3Dmol.createViewer(element);
        // viewer.addSphere({ center: {x:4.76319175, y:-2.95447046, z:15.12600266}, radius: 3.5, color: '#B87333', alpha:0.6 });
        // viewer.addBox({ center: {x:0, y:0, z:0}, dimensions: {w:10, h:10, d:10}, color:'blue', alpha:0.5 });
        viewer.addLine({ start:{x:-10, y:0, z:0}, end:{x:10, y:0, z:0}, color:'red'})
        viewer.addLine({ start:{x:0, y:-10, z:0}, end:{x:0, y:10, z:0}, color:'blue'})
        viewer.addLine({ start:{x:0, y:0, z:-10}, end:{x:0, y:0, z:10}, color:'green'})
        
        viewer.addBox({corner:{x:0,y:0,z:0},dimensions: {w:10,h:10,d:10},color:'magenta',opacity:0.5});
        viewer.addBox({corner:{x:0,y:0,z:0},dimensions: {w:-10,h:-10,d:10},color:'yellow',opacity:0.5});
        
        viewer.addBox({corner:{x:0,y:0,z:0},dimensions: {w:10,h:-10,d:-10},color:'cyan',opacity:0.5});
        viewer.addBox({corner:{x:0,y:0,z:0},dimensions: {w:-10,h:10,d:-10},color:'black',opacity:0.5});
        
        viewer.addBox({corner:{x:0,y:0,z:0},dimensions: {w:10,h:10,d:-10},color:'green',opacity:0.5});
        viewer.addBox({corner:{x:0,y:0,z:0},dimensions: {w:-10,h:10,d:10},color:'blue',opacity:0.5});
        
        viewer.addBox({corner:{x:0,y:0,z:0},dimensions: {w:-10,h:-10,d:-10},color:'red',opacity:0.5});
        viewer.addBox({corner:{x:0,y:0,z:0},dimensions: {w:10,h:-10,d:10},color:'white',opacity:0.5});
        
        viewer.addLabel('X',{position:{x:10, y:0, z:0}})
        viewer.addLabel('Y',{position:{x:0, y:10, z:0}})
        viewer.addLabel('Z',{position:{x:0, y:0, z:10}})

    """

    b = "      viewer.addModel(`" + open("test.xyz").read() + "\n`,'xyz')"

    c = """      
        
            viewer.setStyle({model:-1},{stick:{}})
            viewer.render();
        });
    </script>"""

    with open("vmol.html",'w') as f:
        f.writelines([a,b,c])
        f.close()